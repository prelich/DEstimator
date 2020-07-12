% Sample_Trajectory
% Script to plot a sample trajectory to include in the diffusion estimator
% document
% Author: PKR, UNM, Dec 2015

close all;
clear;
% make this script deterministic
rng(8675309);


%% Create an instance of the Simulation Class and Generate Trajectories
A = DSim; % Call the constructor

% Make 100 random trajectories and pick the most photo-genic one
samples = 1e2;

% figure out the conversion parameters from simulation to experiment
pixelSize = 0.1; % 0.1 um/pixel
frameTime = 0.01; % 0.01 sec/frame  
Dsim = 1e-1; % 1e-1 um^2/sec so that the particle covers a distance

% Determine simulation parameters
% Particle motion and Emission
A.D = Dsim*frameTime/pixelSize^2; % units of px^2/s
A.FrameSize = 100; % maximum frame length for a trajectory
A.MeanPhot = 200; % mean emission rate of photons per frame 
A.PsfSigma = 1; % point spread function sigma in pixels
A.Texp = 1; % exposure time
A.StateParams = [0.2 0.05]; % Kon and Koff states for a `blinking' probe

% Background parameters
A.BGsigma = 5;
A.BGTotalI = 2000;
A.Periodic = 21;

% Number of Trajectories to Sample
A.N = samples;

A.runSim;

% Find the index value with the longest spread for the best image
xdist = A.TrajLimits{1}(3)-A.TrajLimits{1}(1);
ydist = A.TrajLimits{1}(4)-A.TrajLimits{1}(2);
TrajDist = xdist^2 + ydist^2;
index = 1;
for ii = 2:samples
    xdist = A.TrajLimits{ii}(3)-A.TrajLimits{ii}(1);
    ydist = A.TrajLimits{ii}(4)-A.TrajLimits{ii}(2);
     if TrajDist<(xdist^2 + ydist^2)
         index = ii;
         TrajDist = xdist^2 + ydist^2;
     end
end
A.Index = index;
% plot the trajectory
A.plotSim;

%save figures
set(A.Figures.mov,'Units','Inches');
pos = get(A.Figures.mov,'Position');
set(A.Figures.mov,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

set(A.Figures.line,'Units','Inches');
pos = get(A.Figures.line,'Position');
set(A.Figures.line,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(A.Figures.line,'traj_gaussBG_line','fig');
print(A.Figures.line,'traj_gaussBG_line','-dpdf','-r0')
print(A.Figures.line,'traj_GaussBG_line','-dpng');

saveas(A.Figures.mov,'traj_gaussBG_mov','fig');
print(A.Figures.mov,'traj_gaussBG_mov','-dpdf','-r0')
print(A.Figures.mov,'traj_gaussBG_mov','-dpng');
