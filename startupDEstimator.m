%
% Mark J. Olah (mjo@cs.unm.edu)
% 09/2014
%

function startupDEstimator( debug, destimator_path)
    %
    % Sets up the paths to matlab and mex code for DEstimator for windows and linux
    %
    % Inputs:
    %   debug - boolean, determines if we should use debugging version of mex libararies or not (Default: false)
    %   destimator_path - The full path to the DEstimator folder.  Defaults to the location of this setupDEstimator file.
    %

    if nargin==0
        debug=false;
    end
    if nargin<=1
        destimator_path=strsplit(which('startupDEstimator'),'%');% Remove comments from which()
        [destimator_path,~,~]=fileparts(destimator_path{1});
    end

    if ispc()
        if debug
            mex_sub_dir='mex.w64.debug';
        else
            mex_sub_dir='mex.w64';
        end
    elseif isunix()
        if debug
            mex_sub_dir='mex.glnxa64.debug';
        else
            mex_sub_dir='mex.glnxa64';
        end
    else
        error('Platform not supported.');
    end
    destimator_mex_path=fullfile(destimator_path,'mex',mex_sub_dir);
    destimator_matlab_path=fullfile(destimator_path,'matlab');
    addpath(destimator_mex_path);
    addpath(fullfile(destimator_matlab_path, 'utils')); %this is where genpath_safe is
    addpath(genpath_safe(destimator_matlab_path));
end
