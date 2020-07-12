classdef DSim < handle
    %   Class used to generate and simulate diffusing blobs across a non-uniform
    %   background or changing photon count for a trajectory
    %   Development package for Physical Review E
    %   Camera based Diffusion Simulator
    %   Author: Peter K. Relich, UNM November 2015
    
    properties (SetAccess = protected)
        % output values
        X; % Track Coordinates, N+1 positions of [x y t] coords to make N localizations
        O; % Observed Localizations, Times, Photons, BG: N*[o_x, o_y, phot, bg, t]
        CRLB; % CRLB values for the parametrized O; N*[o_x, o_y, phot, bg]
        BinPhot; % Binned photon positions
        PhotStarts; % Boundary values for Binned photons
        BG_Dictionary; % Local Background Dictionary for each of the Binned Photons
        BG_Model; % Padded Background, for display functions
    end
    
    properties
       % general simulation parameters, inputs
       N = 1; % Number of particles to simulate
       FrameSize = 100; % Number of frames per particle
       D=0.1; % Simulated Diffusion Coefficient px^2/frame
       MeanPhot = 200; % expected photons per frame
       PsfSigma=1; % Particle Point Spread Function
       StateParams = [0.2 0.05]; % Kon, Koff rates
       Texp=1; % Exposure Time
    end 
    
    properties
       % background simulation parameters
       Periodic = 21; % Number of pixels before background repeats
       BGsigma = 5; % Sigma of background 'point' source
       BGTotalI = 2000; % Mean bg photons that eminate from a bg 'point' source
    end
    
    properties
       % plotting properties 
       Index = 1; % Value to index into for plotting
       TrajLimits; % positions of photon emissions furthest from center of mass in x,y directions
       ROI = [15 15]; % [x y] dim of pixels to display (centered about a trajectory)
    end
    
    properties (Transient = true, SetAccess = protected)
        % Properties like Emission Times and Intermediate Trajectory Positions
        EmitTime; % The actual emission time of individual photons per track (cell array)
        EmitExp; % Emission exposure time per frame % e.g. a particle that is off will have a time of 0
        AvPos; % Average positions of a particle at a particular emission time
        BG_Val; % Model calculated background value of AvPos locations (prior to noise on BG)
        BG_Stack; % stack of background images generated for plot function
        Data; % Output data set
        Figures; % Figure plot handles and axes
    end
    
    properties (SetAccess = protected)
        States; % Continuous Emission State Transitions per track (cell array)
    end
    
    methods
        % Constructor
        % Very simple, it doesn't need input arguments.
        function obj = DSim(varargin)

        end
        
        % parameter restrictions
        function set.N(obj, value)
           if (value < 1 || mod(value,1))
               error('Minimum N particles is 1 and must be an integer value');
           else
               obj.N = value;
           end
        end
        
        function set.Texp(obj, value)
           if (value < 0 || value > 1)
               error('Exposure time must be between 0 and 1');
           else
               obj.Texp = value;
           end
        end
        
        function set.StateParams(obj, value)
           if (any(value < 0) || any(value > 1))
               error('Emission rates must be between 0 and 1');
           else
               obj.StateParams = value;
           end
        end
        
        function set.MeanPhot(obj, value)
            if value < 0
                error('Mean photon rates must be postive')
            else
                obj.MeanPhot = value;
            end
        end
        
        % Main Simulation Running Script
        function runSim(obj)
           obj.genTrueDiffusion; % create pure diffusion trajectories
           obj.expStates; % determine emission states
           obj.genEmitTimes; % get emission times during periods of on states
           obj.getPhotonPos; % get the binned photons and average position of pure diffusion process
           obj.gaussBG; % build the gaussian lattice background
           obj.returnBGval; % returns a realization of the fit back ground value
           obj.getCRLB; % returns the expected uncertainties of the MLE vectors
           obj.getObs; % generates the MLE vectors for the fit positions, photons, and background with CRLB
        end
        
        
        % Generate trajectories, put them in cell arrays
        function genTrueDiffusion(obj)
            starts = obj.Periodic*rand(2,obj.N);
            starts = reshape(starts,[1 2 obj.N]);
            Xtemp = sqrt(2*obj.D)*cumsum(randn(obj.FrameSize+1,2,obj.N));
            Xtemp = Xtemp + repmat(starts,obj.FrameSize+1,1,1);
            if obj.N > 1
                frameAssoc = (1:obj.FrameSize+1)';
                frameAssoc = repmat(frameAssoc,[1 1 obj.N]);
                obj.X = mat2cell([Xtemp frameAssoc],obj.FrameSize+1,3,ones(1,obj.N));
                obj.X = reshape(obj.X,[],1);
            else
                obj.X = mat2cell([Xtemp (1:obj.FrameSize+1)'],obj.FrameSize+1,3)';
            end
        end
                
        
        function getCRLB(obj)
            % written as a vector operation
            % for a given set of Avpos, photon and BG counts, return the CRLB
            obj.CRLB = cell(obj.N,1);
            for ii = 1:obj.N
                % build the theta vector
                % shift the x and y values
                shift = [obj.AvPos{ii}(:,1),obj.AvPos{ii}(:,2)];
                sz = 6*obj.PsfSigma+1;
                shift = round(shift-sz/2);
                box_coord = obj.AvPos{ii}(:,1:2)-shift;

                tempEmitExp = obj.EmitExp{ii};
                tempEmitExp = tempEmitExp(obj.AvPos{ii}(:,4)); % only take frames that generated photons!
                theta = [box_coord, obj.MeanPhot*tempEmitExp, obj.BG_Val{ii}]';
                Nfits = size(theta,2);
                obj.CRLB{ii} = obj.gaussianCRLB(theta,obj.PsfSigma,round(sz),Nfits);
            end
            
        end
        
        function [Obs, Var] = filterCRLB(obj)
            % function that removes all observations with CRLB values that
            % are over a half pixel (what is typically thrown out for localizations)
            Obs = obj.O;
            Var = obj.CRLB;
            for ii = 1:obj.N
               remind = (obj.CRLB{ii}(:,1) > 0.5 | obj.CRLB{ii}(:,2) > 0.5);
               Var{ii}(remind,:) = [];
               Obs{ii}(remind,:) = [];
            end
        end
        
        function getObs(obj)
            % function to return the Obs vectors, used to find CRLB
            obj.O = cell(obj.N,1);
            for ii = 1:obj.N
                Ox = obj.AvPos{ii}(:,1)+sqrt(obj.CRLB{ii}(:,1)).*randn(length(obj.CRLB{ii}(:,1)),1);
                Oy = obj.AvPos{ii}(:,2)+sqrt(obj.CRLB{ii}(:,2)).*randn(length(obj.CRLB{ii}(:,2)),1);
                obj.O{ii} = [Ox Oy obj.AvPos{ii}(:,3) ...
                    obj.BG_Val{ii} obj.AvPos{ii}(:,4)];
            end
        end
        
        
        function returnBGval(obj)
            % returns the estimated background from averaging a 3*psfsigma
            % + 1 region about some point in the ROI, Associates averaged
            % positions with the surrounding BG dictionary as well as the
            % corresponding poisson realization.
            obj.BG_Val = cell(obj.N,1);
            for ii = 1:obj.N
                BG_pixel = round(obj.AvPos{ii}(:,1:2)+1/2);
                BG_loc = mod(BG_pixel-1,obj.Periodic)+1;
                obj.BG_Val{ii} = zeros(length(obj.AvPos{ii}(:,1)),1);
                for jj = 1:length(obj.AvPos{ii}(:,1))
                    obj.BG_Val{ii}(jj) = obj.BG_Dictionary(BG_loc(jj,1),BG_loc(jj,2));
                end
            end           
        end
        
        % Generate emission states
        function expStates(obj)
            % generates on and off states for each probe per frame given
            % an on and off rate parameterized by exponential distributions.
            % a function for simulating organic dyes.
            T = obj.FrameSize;
            % estimate 5 std of transitions from on and off rates
            Kon = obj.StateParams(1);
            Koff = obj.StateParams(2);
            % back of the envelope calculation of the inverse of average life time per
            % state * frames - 1
            avrate = 1/(1/Kon + 1/Koff)*2*(T-1);
            % there will be 2x samples, partitioned into halves per state
            samples = ceil((avrate+1 + 5*sqrt(avrate+1)/2));
            Pon = Kon/(Kon+Koff); % get a steady state start probability
            StateStart = rand(obj.N,1)<Pon; % 1 is the on state, 0 is the off state
            % generate State Transition Cell
            StateT = cell(obj.N,1);
            % generate emission exposure time per frame cell
            obj.EmitExp = cell(obj.N,1);
            K = [Kon Koff]; % put the states into a vector
            for ii = 1:obj.N
                % generate random numbers
                CDFval = rand(2*samples-1,1);
                % generate a K switch vector based on the start state
                Kswitch = [ones(1,samples)*K(2-StateStart(ii)); ones(1,samples)*K(StateStart(ii)+1)];
                Kswitch = Kswitch(:);
                % get a list of the states per transition
                Kstates = [ones(1,samples)*StateStart(ii); ones(1,samples)*mod(StateStart(ii)+1,2)];
                Kstates = Kstates(:);
                % convert the CDFvalues to the wait times and do a
                % cumulative sum to get state switch times
                xx = [0 ; -log(1-CDFval)./Kswitch(2:end)];
                xx = cumsum(xx); % get change times
                % remove elements longer than the simulation length
                remind = xx > T;
                xx(remind) = [];
                Kstates(remind) = [];
                Kswitch(remind) = [];
                % generate the state matrix of transition times and states
                StateT{ii} = [xx Kstates Kswitch];
                
                % get the exposure time per frame given an on emission state
                k_index = mod(StateT{ii}(end,2),2)+1;
                % the last state is always off!
                tempStates = [StateT{ii}; obj.FrameSize, 0, obj.StateParams(k_index)];
                exp_frame = zeros(T,1);
                counter = 0;
                while counter < size(StateT{ii},1)
                    counter = counter+1;
                    state_p = tempStates(counter,:);
                    state_n = tempStates(counter+1,:);
                    % determine if the exposure time is 0, 1, or fractional!
                    exp_frame(ceil(state_p(1)+1):floor(state_n(1))) = ...
                        exp_frame(ceil(state_p(1)+1):floor(state_n(1))) + state_p(2);
                    fractime = min(mod(state_n(1),1),state_n(1)-state_p(1));
                    exp_frame(ceil(state_n(1))) = exp_frame(ceil(state_n(1))) ...
                        + fractime*state_p(2) + (1-fractime)*state_n(2);
                end
                obj.EmitExp{ii} = exp_frame;
            end
            obj.States = StateT;
        end
        
        function gaussBG(obj)
            
            % Draw the Gaussian model
            obj.BG_Model = DSim.finitegausspsf(obj.Periodic,obj.BGsigma,obj.BGTotalI,0,[obj.Periodic/2 obj.Periodic/2]);
            
            NeighborSearch = round(3*obj.PsfSigma);
            % Pad the BGModel with the neighbor search region to get averaged
            % background values per central pixel location of the gaussian emitter
            BG_Padded = zeros(size(obj.BG_Model)+2*NeighborSearch);
            BG_Padded(NeighborSearch+(1:size(obj.BG_Model,1)),NeighborSearch+(1:size(obj.BG_Model,2))) = obj.BG_Model;
            % Adjacent Neighbors
            BG_Padded(1:NeighborSearch,NeighborSearch+(1:size(obj.BG_Model,2))) = ...
                obj.BG_Model((end-NeighborSearch+1):end,:);
            BG_Padded((end-NeighborSearch+1):end,NeighborSearch+(1:size(obj.BG_Model,2))) = ...
                obj.BG_Model(1:NeighborSearch,:);
            BG_Padded(NeighborSearch+(1:size(obj.BG_Model,1)),1:NeighborSearch) = ...
                obj.BG_Model(:,(end-NeighborSearch+1):end);
            BG_Padded(NeighborSearch+(1:size(obj.BG_Model,1)),(end-NeighborSearch+1):end) = ...
                obj.BG_Model(:,1:NeighborSearch);
            % Diagonal Neighbors
            BG_Padded(1:NeighborSearch,1:NeighborSearch) = ...
                obj.BG_Model((end-NeighborSearch+1):end,(end-NeighborSearch+1):end);
            BG_Padded(1:NeighborSearch,(end-NeighborSearch+1):end) = ...
                obj.BG_Model((end-NeighborSearch+1):end,1:NeighborSearch);
            BG_Padded((end-NeighborSearch+1):end,1:NeighborSearch) = ...
                obj.BG_Model(1:NeighborSearch,(end-NeighborSearch+1):end);
            BG_Padded((end-NeighborSearch+1):end,(end-NeighborSearch+1):end) = ...
                obj.BG_Model(1:NeighborSearch,1:NeighborSearch);
            % Figure out the dictionary Background given the search region of pixels
            obj.BG_Dictionary = zeros(size(obj.BG_Model));
            % create dictionary of averaged background values over search region
            SchRgn = 2*NeighborSearch;
            for ii = 1:size(obj.BG_Model,1)
                for jj = 1:size(obj.BG_Model,2)
                    obj.BG_Dictionary(ii,jj) = mean(mean(BG_Padded(ii:(ii+SchRgn),jj:(jj+SchRgn))));
                end
            end
            
            
        end
        
        % Generate Emission Times
        function genEmitTimes(obj)
            % function to generate photon emission times continuously
            % during a particle's trajectory, follows a poisson process
            obj.EmitTime = cell(obj.N,1);
            % loop over each frame and determine photon emission times
            for ii = 1:obj.N
                tempEmit = cell(obj.FrameSize);
                for jj = 1:obj.FrameSize % loop frame times
                    % sample out 5 sigmas of random numbers from the max rate
                    % speed is more important than memory management here
                    sample = ceil(obj.MeanPhot+5*sqrt(obj.MeanPhot));
                    CDFval = rand(1,sample); % get all the emission CDF times
                    % convert the CDFvalues to the wait times and do a
                    % cumulative sum to get emission times
                    xx = -log(1-CDFval)/obj.MeanPhot;
                    xx = cumsum(xx);
                    % remove anything greater than the exposure time
                    xx(xx > obj.Texp) = [];
                    tempEmit{jj} = xx+jj-1; % frame 1 starts at time 0
                end
                obj.EmitTime{ii} = [tempEmit{:}]';
            end
            % remove emission times that occur during an off state
            if ~isempty(obj.States)
                for ii = 1:obj.N % loop over probed tracks
                    % put in a 0 state at the end of the trajectory, always off!
                    k_index = mod(obj.States{ii}(end,2),2)+1;
                    tempStates = [obj.States{ii}; obj.FrameSize, 0, obj.StateParams(k_index)];
                    
                    transitions = size(obj.States{ii},1);
                    remind = zeros(length(obj.EmitTime{ii}),1);
                    for jj = 1:transitions
                        % remove emissions that occur during an off state
                        if tempStates(jj,2) == 0
                            remind = remind | (obj.EmitTime{ii} > tempStates(jj,1) ...
                                & obj.EmitTime{ii} < tempStates(jj+1,1));
                        end
                    end
                    if nnz(remind)
                        obj.EmitTime{ii}(remind) = [];
                    end
                end
            end
        end

        % function to calculate photon positions from emission times and
        % coordinates with the brownian bridge formalism.
        % photons are binned to pixels and stored in aggregate
        function getPhotonPos(obj)
           % 2-D cell array of sparse matrices required to draw sub-images
           obj.BinPhot = cell(obj.FrameSize,obj.N);
           % Absolute values to start positions required to draw global-images
           obj.PhotStarts = cell(obj.FrameSize,obj.N);
           % Averaged position of the counted photons (equivalent to MLE)
           obj.AvPos = cell(obj.N,1);
           % Determine the span of a trajectory
           obj.TrajLimits = cell(obj.N,1);
           % loop over each emitter and determine photon positions
           for ii = 1:obj.N
               count = 0; % start the counter for intermittency
               timevec = obj.EmitTime{ii};
               tempAvPos = zeros(obj.FrameSize,4);
               obj.TrajLimits{ii} = inf(4,1); % min at inf
               obj.TrajLimits{ii}(3:4) = -inf; % max at -inf
               for jj = 1:obj.FrameSize % loop over frames
                   timeInterval = timevec(timevec>jj-1 & timevec<jj)-jj+1;
                   PhotCounts = length(timeInterval);
                   if isempty(timeInterval) % skip intermittent points
                       continue;
                   end
                   count = count+1;
                   frameInput = [obj.X{ii}(jj,:) obj.X{ii}(jj+1,:)];
                   % calculate the positions due to an interpolated
                   % brownian bridge
                   intpos = bBridge(obj,timeInterval,frameInput);                   
                   % add diffraction based randomness
                   xx = intpos(:,1)+obj.PsfSigma*randn(PhotCounts,1);
                   yy = intpos(:,2)+obj.PsfSigma*randn(PhotCounts,1);
                   % Figure out end points of photon emission for a trajectory
                   minLimit = min(obj.TrajLimits{ii}(1:2),[min(xx);min(yy)]);
                   maxLimit = max(obj.TrajLimits{ii}(3:4),[max(xx);max(yy)]);
                   obj.TrajLimits{ii} = [minLimit;maxLimit];
                   % bin photons to their respective pixels
                   % store as a sub-matrix with start x and start y
                   tempbin = round([xx yy]);
                   startx = min(tempbin(:,1));
                   starty = min(tempbin(:,2));
                   tempbin = [tempbin(:,1)-startx+1 tempbin(:,2)-starty+1];
                   photmat = sparse(tempbin(:,1),tempbin(:,2),1);
                   % store the sparse photon matrix, as well as the starts
                   obj.BinPhot{jj,ii} = photmat;
                   obj.PhotStarts{jj,ii} = [startx,starty];
                   % store the averaged positions for CRLB and Diffusion
                   % Calculations, average is of binned pixels
                   tempAvPos(count,:) = [mean(intpos,1) PhotCounts jj]; % store frame location
               end
               tempAvPos(count+1:end,:) = []; % delete unused vector points
               % Store the averaged positions
               obj.AvPos{ii} = tempAvPos;
           end            
        end
        
        % function to generate the interpolated positions with the brownian
        % bridge formalism
        function intpos = bBridge(obj,timevec,frameInput)
            % delt is the time gap of the bridge, A and B are the bridge
            % end points
            delt = frameInput(6)-frameInput(3);
            A = [frameInput(1) frameInput(2)];
            B = [frameInput(4) frameInput(5)];
            % intermediate variables to speed up bridge calculations
            dL = (diff([0; timevec]));
            dR = delt-timevec;
            dN = delt-[0; timevec(1:end-1)];
            alpha = dL./dN;
            % bridge covariances in X and Y dimensions
            bvar(:,1) = sqrt(2*obj.D*alpha.*dR).*randn(length(alpha),1);
            bvar(:,2) = sqrt(2*obj.D*alpha.*dR).*randn(length(alpha),1);
            % generate interpolated positions
            intpos = zeros(length(timevec),2);
            intpos(1,:) = A*(1-alpha(1)) + alpha(1)*B + bvar(1,:);
            for ii = 2:length(timevec)
                intpos(ii,:) = intpos(ii-1,:)*(1-alpha(ii)) + alpha(ii)*B + bvar(ii,:);
            end
        end
        
        
        % Function to plot the trajectory of a given index value over the
        % gaussian lattice
        function plotSim(obj)
            obj.BG_Stack = zeros([obj.ROI obj.FrameSize]);
            modelBG = zeros([obj.ROI]);
            % determine what the ROI corresponds to in absolute position
            div = [obj.ROI(1)/2, obj.ROI(2)/2];
            % find the center of the trajectory limits and ROI edges
            trajcenter = [(obj.TrajLimits{obj.Index}(1)+obj.TrajLimits{obj.Index}(3))/2 ...
                (obj.TrajLimits{obj.Index}(2)+obj.TrajLimits{obj.Index}(4))/2];
            ROIlim = round(trajcenter - div);
            
            % loop over all BG pixels and associate the proper BG model value
            for ii = 1:obj.ROI(1)
                for jj = 1:obj.ROI(2)
                    ROItoBG_X = mod(ROIlim(1)+ii-1,20)+1;
                    ROItoBG_y = mod(ROIlim(2)+jj-1,20)+1;
                    modelBG(ii,jj) = obj.BG_Model(ROItoBG_X,ROItoBG_y);
                end
            end
            
            % loop over each frame to perform poisson noise realizations on the bg stack
            for ii = 1:obj.FrameSize
                obj.BG_Stack(:,:,ii) = poissrnd(modelBG);
            end
            
            % loop over each frame and added photons to the BG stack
            % tempIm = zeros(obj.ROI);
            for ii = 1:obj.FrameSize
                % convert Binned Photons to Global Pixels
                obj.Data(:,:,ii) = obj.BG_Stack(:,:,ii);
                % Figure out where the photons should fit globally
                PhotEnds = obj.PhotStarts{ii,obj.Index};
                if isempty(PhotEnds)
                    continue;
                end
                % Paint the particle according to its relative position in the ROI
                [iphot, jphot, vphot] = find(obj.BinPhot{ii,obj.Index});
                iphot = iphot + PhotEnds(:,1) - ROIlim(1);
                jphot = jphot + PhotEnds(:,2) - ROIlim(2);
                for kk = 1:length(iphot)
                    if iphot(kk)>obj.ROI(1) || jphot(kk)>obj.ROI(2) || iphot(kk)<1 || jphot(kk)<1
                        continue; % can't draw outside of the ROI
                    end
                    obj.Data(iphot(kk),jphot(kk),ii) = obj.Data(iphot(kk),jphot(kk),ii) + vphot(kk);
                end
            end
            
            % plot something!
            generateFigPlots(obj,ROIlim);
            
        end
        
        
        
        % plotting function for this method
        function generateFigPlots(obj,ROIlim)
            
            % first plot the movie
            dipshow(obj.Data);
            hold on;
            % only plots a single trajectory over a corresponding bg region
            val = obj.Index;
            plot(obj.O{val}(:,2)-ROIlim(2),obj.O{val}(:,1)-ROIlim(1),'Color',[0.9 0.7 0],'LineWidth',2);
            obj.Figures.mov = gcf;
            obj.Figures.movha = gca;
            
            plot([1; 6], [3; 3], '-w',  'LineWidth', 2);
            movh2 = text(5,4, '500 nm', 'HorizontalAlignment','right');
            movh2.Color = [1 1 1];
            movh2.FontWeight = 'bold';
            movh2.Parent.Parent.Color = [1 1 1];
            
            hold off;
            % then plot the histogram of the first trajetory
            figure;
            pixelCRLB = obj.CRLB{val}(:,1);
            plot(obj.O{val}(:,5),sqrt(pixelCRLB(:,1)),'o');
            obj.Figures.line = gcf;
            obj.Figures.lineha = gca;
            obj.Figures.lineha.FontSize = 12;
            obj.Figures.lineha.FontWeight = 'bold';
            % obj.Figures.lineha.Title.String = 'Localization \sigma at Each Frame';
            obj.Figures.lineha.XLabel.String = 'Frame Index';
            obj.Figures.lineha.YLabel.String = 'Estimated Localization \sigma (px)';
            obj.Figures.lineha.Children.LineWidth = 3;
            obj.Figures.lineha.Children.Color = [0 0 0];
        end
        
    end
    
    
    % Methods borrowed from Lidke Lab repository
    methods (Static = true)
        function [out] = finitegausspsf(Npixels,sigma,I,bg,cor)
            %finiteGaussPSFerf   Make Gaussian Spots using finite pixel size
            %
            %   [out] = finiteGaussPSFerf(Npixels,sigma,I,bg,cor)
            %       INPUT
            %   Npixels:    linear size in pixels
            %   sigma:      PSF sigma in pixels, scaler gives symmetric, [sx sy] gives
            %               asymmetric.
            %   I:          Photons/frame
            %   bg:         background/pixel
            %   cor:        coordinates array size of [2 N] where N is number of spots
            %               to generate.
            %       OUTPUT
            %   out:        3D stack of images.
            if (nargin < 5)
                cor = [(Npixels-1)/2 (Npixels-1)/2 0];
            end
            if (nargin <4)
                error('Minimal usage: finitegausspsf(Npixels,sigma,I,bg)');
            end            
            if (bg == 0)
                bg = 10^-10;
            end
            Ncor=size(cor,1);
            
            x=repmat((0:Npixels-1)',[1 Npixels]);
            y=x';
            X=repmat(x,[1 1 Ncor]);
            Y=repmat(y,[1 1 Ncor]);
            
            Xpos=repmat(shiftdim(cor(:,1),-2),[Npixels,Npixels,1]);
            Ypos=repmat(shiftdim(cor(:,2),-2),[Npixels,Npixels,1]);
            
            if size(I,1) > 1
                I=repmat(shiftdim(I,-2),[Npixels,Npixels,1]);
                if max(size(Xpos) ~= size(I))
                    error('Size of I and maybe others are incorrect.');
                end
            end
            if size(bg,1) > 1
                bg=repmat(shiftdim(bg,-2),[Npixels,Npixels,1]);
                if max(size(Xpos) ~= size(bg))
                    error('Size of bg and maybe others are incorrect.');
                end
            end
            if length(sigma)==1
                sigmay = sigma;
                sigmax = sigma;
            else
                sigmax = sigma(1);
                sigmay = sigma(2);
            end
            gausint=I/4.*((erf((X-Xpos+.5)./(sqrt(2)*sigmax))-erf((X-Xpos-.5)./(sqrt(2)*sigmax))).*...
                (erf((Y-Ypos+.5)./(sqrt(2)*sigmay))-erf((Y-Ypos-.5)./(sqrt(2)*sigmay))))+bg;
            out=gausint;
        end

        % function to estimate the theoretical Fisher Positions given "known"
        % probe information (photons, background)
        function CRLB = gaussianCRLB(theta,PSFSigma,sz,Nfits)
            % Initializee fisher matrix and CRLB data stacks
            CRLB = zeros(Nfits,4); % standard 4 parameter fitting
            M = zeros(4,4,Nfits); 
            Minv = zeros(4,4,Nfits);            
            % Calculating the CRLB
            for ii = 1:sz
                for jj = 1:sz
                    PSFx = DSim.IntGauss1D( ii-1, theta(1,:), PSFSigma );
                    PSFy = DSim.IntGauss1D( jj-1, theta(2,:), PSFSigma );                    
                    model = theta(4,:) + theta(3,:).*PSFx.*PSFy;
                    % Calculating Derivatives
                    dudt(1,:) = DSim.DerivativeIntGauss1D( ii-1, theta(1,:), PSFSigma, theta(3,:), PSFy );
                    dudt(2,:) = DSim.DerivativeIntGauss1D( jj-1, theta(2,:), PSFSigma, theta(3,:), PSFx );
                    dudt(3,:) = PSFx.*PSFy;
                    dudt(4,:) = ones(1,Nfits);
                    % Building the Fisher Information Matrix
                    for nn = 1:Nfits
                        M(:,:,nn) = M(:,:,nn) + dudt(:,nn) * dudt(:,nn)' / model(nn);
                    end
                end
            end
            % Matrix inverse (CRLB = F^-1) and output assignments
            for nn = 1:Nfits
                Minv(:,:,nn) = inv(M(:,:,nn));
            end
            % Write to global arrays
            for nn = 1:Nfits
                for kk = 1:4
                    CRLB(nn,kk) = Minv(kk,kk,nn);
                end
            end
        end
        
        function [dudt, d2udt2] = DerivativeIntGauss1D( ii, x, sigma, N, PSFy )
            % \brief computes the 1st and 2nd derivatives of the 1D gaussian
            %
            % Inputs:
            % \ii index location (used as an off-set) from the external loop
            % \x x-coordinate
            % \sigma sigma
            % \N number of Photons
            % \PSFy Point Spread Function coefficient that remains constant in this derivative
            %
            % Outputs:
            % \Idudt first derivative
            % \Id2udt2 second derivative            
            a = exp(-0.5*( (ii + 0.5 - x)./sigma ).^2);
            b = exp(-0.5*( (ii - 0.5 - x)./sigma ).^2);
            dudt = -N.*PSFy.*(a - b) ./ (sqrt(2.0*pi).*sigma);
            d2udt2 = -N.*( (ii + 0.5 - x).*a - (ii - 0.5 - x).*b ).*PSFy ./ ( sqrt(2.0*pi).*sigma.^3 );
        end
        
        function AvInt = IntGauss1D( ii, x, sigma )
            % \brief computes average of off-centered error function integrals
            %
            % Inputs:
            % \ii counting parameter, sets off-set in error function
            % \x position parameter in error function
            % \sigma sigma value of the PSF
            %
            % Output: avInt, the average error function integral value
            norm = 0.5 ./ (sigma.^2);
            AvInt = 0.5 * ( erf( (ii - x + 0.5).*sqrt(norm) ) - erf( (ii - x -0.5).*sqrt(norm) ) );
        end
    end
    
end


