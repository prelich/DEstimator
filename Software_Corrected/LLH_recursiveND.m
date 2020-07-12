function LLH = LLH_recursiveND(D, Obs, T, SE, exposureT)
% This is an extra function to generate the log likelihood of a
% single D value given a multi-dimensional trajectory presented in a matrix format.
%
% Usage: llH = LLH_recursiveND(D, Obs, T, SE,
% exposureT;
%
%   INPUTS:
%       D - Vector: sampled diffusion constant
%       Obs - size=Nxdim: observed points
%       T - size=N: time vector (must be incremental)
%       SE - size=Nxdim: observed localization errors
%       exposureT - the exposure time.  REQUIRED.
LLH = zeros(length(D),1);
% loop over each trajectory to get the log likelihood value
for ii = 1:size(Obs,2)
    tObs = Obs(:,ii);
    tSE = SE(:,ii);
    llh_temp = LLH_recursive1D(D, tObs, T, tSE, exposureT);
    LLH = LLH + llh_temp;
end
end