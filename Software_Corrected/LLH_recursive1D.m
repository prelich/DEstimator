function LLH = LLH_recursive1D(D, Obs, T, SE, exposureT)
% Recursive solution, Matlab implementation
%   INPUTS:
%       D - size=ND: sample diffusion values vector or scalar
%       Obs - size=N: observed points
%       T - size=N: time vector (must be incremental)
%       SE - size=N: observed localization errors
%       exposureT - size=1: Exposure Time. REQUIRED
%   OUTPUT:
%       LLH - size=ND: log likelihood of observing Obs with standard error SE
%       given true diffusion constant is D.
N = length(Obs);
ND = length(D);
dt = diff(T);
LLH = zeros(ND,1);
for nn = 1:ND
    vD = 2*D(nn)*dt;
    vM = SE.*SE - 2*D(nn)*exposureT/6;
    eta = vM(1)+vD(1);
    mu = Obs(1);
    for k = 2:N-1
        alpha = eta+vM(k);
        LLH(nn) = LLH(nn) + log(alpha)+ (Obs(k)-mu)^2/alpha;
        mu = (mu*vM(k)+(Obs(k))*eta)/alpha;
        eta = vD(k)+vM(k)*eta/alpha;
    end
    alpha = eta+vM(N);
    LLH(nn) = LLH(nn) + (N-1)*log(2*pi) + log(alpha) + (Obs(N)-mu)^2/alpha;
end
LLH = -0.5*LLH;  %factored out -0.5 from everything above
end