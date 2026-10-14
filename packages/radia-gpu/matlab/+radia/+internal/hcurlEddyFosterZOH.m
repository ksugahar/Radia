function [adDiagonal, bdGain] = hcurlEddyFosterZOH(decayRates, sampleTime_s)
%HCURLEDDYFOSTERZOH Exact zero-order hold of the diagonal Foster system.
%   For z_dot = -lambda*z + b*u the hold gives z+ = exp(-lambda*T)*z +
%   g*b*u with g = (1 - exp(-lambda*T))/lambda, and g = T for lambda = 0.

arguments
    decayRates double {mustBeFinite, mustBeNonnegative}
    sampleTime_s (1,1) double {mustBeFinite, mustBePositive}
end

lambda = double(decayRates(:));
adDiagonal = exp(-lambda * sampleTime_s);
bdGain = repmat(sampleTime_s, size(lambda));
nonzero = lambda > 0;
bdGain(nonzero) = -expm1(-lambda(nonzero) * sampleTime_s) ./ lambda(nonzero);
end
