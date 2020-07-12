function Information = recursiveFisher1D(D, T, SE, exposureT)
% Recursive solution to the Fisher Information,  Matlab implementation,
% Returns a vector of Fisher information values for a given D, Recursive solution
%   INPUTS:
%       D - size=ND: sample MLE values vector or scalar
%       T - size=N: time vector (must be incremental)
%       SE - size=N: observed localization errors
%       exposureT - size=1: Exposure Time. REQUIRED
%   OUTPUT:
%       Information - size=ND: Fisher information sampled at each potential
%       MLE (D), given the conditional vectors SE and T and sclar exposureT
N = length(SE);
ND = length(D);
dt = diff(T);
Information = zeros(ND,1);
for nn = 1:ND
    vD = 2*D(nn)*dt;
    vM = SE.*SE - 2*D(nn)*exposureT/6;
    % initial expectation values
    expdmu2 = 0;
    % initial derivative values
    dvD = 2*dt;
    dvM = -exposureT/3*ones(length(SE),1);
    dalpha = dvM(1)+dvD(1)+dvM(2);
    % initial recursive values
    alpha = vM(1)+vD(1)+vM(2);
    for k = 2:N-1
        Information(nn) = Information(nn) + (dalpha^2)/2/alpha^2 ...
            + expdmu2/alpha;
        % propagate expectation values, higher values first
        expdmu2 = 1/alpha*(dalpha*vM(k)/alpha - dvM(k))^2 ...
            + (vM(k)/alpha)^2*expdmu2;
        % propagate recursive derivatives, higher derivatives first.
        dalpha = dvD(k)+dvM(k+1)+dvM(k)-2*vM(k)*dvM(k)/alpha ...
            +(vM(k)/alpha)^2*dalpha;
        alpha = vD(k)+vM(k+1)+vM(k)-vM(k)^2/alpha;
    end
    Information(nn) = Information(nn) + (dalpha^2)/2/alpha^2 ...
        + expdmu2/alpha;
end
end