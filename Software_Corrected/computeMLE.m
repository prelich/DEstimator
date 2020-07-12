function [MLE,MLE_LLH]=computeMLE(Obs, T, SE, exposureT)
% Author: MJO, UNM, 2015
% BRIEF: estimates Diffusion constant (isotropic) from a multi-dimensional trajectory
%
% Usage:
%   [MLE,MLE_LLH]=computeMLE(D,Obs, T, SE, exposureT);
%
%   INPUTS:
%       Obs - size=Nxdim: observed points
%       T - size=N: time vector (must be incremental)
%       SE - size=Nxdim: observed localization errors
%       exposureT - size=1: Exposure Time. REQUIRED
%   OUTPUT:
%       LLH - size=ND: log likelihood of observing Obs with standard error SE
%       given true diffusion constant is D.
%
minD = 1e-8;
maxD = 1e8;
maxEvals = 1000;
opts = optimset('MaxFunEvals', maxEvals, 'MaxIter', maxEvals);
[MLE, MLE_LLH] = fminbnd(@(d) -LLH_recursiveND(d, Obs, T, SE, exposureT),minD, maxD,opts);
MLE_LLH = -MLE_LLH; %Change from minimization to maximization
end