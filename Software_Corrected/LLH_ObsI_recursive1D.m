function [LLH, ObsInformation] = LLH_ObsI_recursive1D(D, Obs, T, SE, exposureT)
% Recursive solution, pure Matlab implementation, returns information
% vector as well, only for 1-D trajectories.  Information is additive along
% dimensions (can be called multiple times and summed for each dimension).
%   INPUTS:
%       D - size=ND: sample diffusion values vector or scalar
%       Obs - size=N: observed points
%       T - size=N: time vector (must be incremental)
%       SE - size=N: observed localization errors
%       exposureT - size=1: Exposure Time. REQUIRED
%   OUTPUT:
%       LLH - size=ND: log likelihood of observing Obs with standard error SE
%       given true diffusion constant is D.
%       ObsInformation - size=ND: observed information returned at each D
%       sampled.  Is only guaranteed to be positive and meaningful at the MLE
N = length(Obs);
ND = length(D);
dt = diff(T);
LLH = zeros(ND,1);
ObsInformation = zeros(ND,1);
for nn = 1:ND
    vD = 2*D(nn)*dt;
    vM = SE.*SE - 2*D(nn)*exposureT/6;
    % initial derivative values
    dvD = 2*dt;
    dvM = -exposureT/3*ones(length(SE),1);
    ddmu = 0;
    ddalpha = 0;
    dalpha = dvM(1) + dvD(1) + dvM(2);
    dmu = 0;
    % initial recursive values
    eta = vM(1)+vD(1);
    mu = Obs(1);
    for k = 2:N-1
        alpha = eta+vM(k);
        LLH(nn) = LLH(nn) + log(alpha)+ (Obs(k)-mu)^2/alpha;
        ObsInformation(nn) = ObsInformation(nn) + (ddalpha-dalpha^2/alpha)/alpha/2 ...
            + (Obs(k)-mu)^2/alpha^2*( (dalpha)^2/alpha - ddalpha/2) + ...
            (Obs(k)-mu)/alpha*(2*dalpha*dmu/alpha - ddmu) + (dmu)^2/alpha;
        % propagate recursive derivatives, higher derivatives first.
        ddmu = ddmu*vM(k)/alpha + 2/alpha*dmu*(dvM(k) - vM(k)*dalpha/alpha)...
            + (Obs(k)-mu)/alpha^2*(2*dvM(k)*dalpha - 2*vM(k)*dalpha^2/alpha ...
            + ddalpha*vM(k));
        ddalpha = (vM(k)/alpha)^2*ddalpha + 4*vM(k)/alpha^2*dvM(k)*dalpha ...
            - 2*(vM(k)/alpha)^2/alpha*dalpha^2 - 2/alpha*(dvM(k))^2;
        dmu = dmu*vM(k)/alpha+ (Obs(k)-mu)/alpha*(vM(k)/alpha*dalpha-dvM(k));
        dalpha = dvD(k) + dvM(k+1) + dvM(k) + vM(k)/alpha*(vM(k)/alpha*dalpha - 2*dvM(k));
        % update recursive variables
        mu = (mu*vM(k)+(Obs(k))*eta)/alpha;
        eta = vD(k)+vM(k)*eta/alpha;        
    end
    alpha = eta+vM(N);
    LLH(nn) = LLH(nn) + (N-1)*log(2*pi) + log(alpha) + (Obs(N)-mu)^2/alpha;
    ObsInformation(nn) = ObsInformation(nn) + (ddalpha-dalpha^2/alpha)/alpha/2 ...
            + (Obs(N)-mu)^2/alpha^2*( (dalpha)^2/alpha - ddalpha/2) + ...
            (Obs(N)-mu)/alpha*(2*dalpha*dmu/alpha - ddmu) + (dmu)^2/alpha;
end
LLH = -0.5*LLH;  %factored out -0.5 from everything above
end