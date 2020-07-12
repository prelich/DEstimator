function Deviation = waldInterval( ObsInf, pval )
% WaldInterval returns the Wald Interval for a given MLE and Observed
% Information
%   The script is basis agnostic, but the Wald Interval itself is not
%   invariant to transformation, so make sure that the Observed
%   Information is in the correct basis or the interval will be meaningless
%
%   Inputs: ObsInf - Observed Information at the MLE, must be positive.
%           pval - % of the Gaussian to exclude from the interval
%
%   Outputs: Deviation - Deviation from the MLE in +/- of the corresponding
%   basis to establish the wald interval
Deviation = sqrt(2/ObsInf)*erfinv(1-pval);

end

