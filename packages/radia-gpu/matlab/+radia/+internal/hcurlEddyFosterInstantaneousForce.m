function force_N = hcurlEddyFosterInstantaneousForce(forceOperator, modalState, coilCurrent)
%HCURLEDDYFOSTERINSTANTANEOUSFORCE Instantaneous reduced Lorentz force.
%   For real time-domain signals the force is F_k(t) = sum K(k,j,b)*z_j(t)*i_b(t).
%   Its cycle average equals the phasor time average
%   0.5*real(sum K*Z*conj(I)) returned by evaluateHCurlEddyFosterForce, so
%   the time-domain block must not apply the phasor factor 0.5.

nState = size(forceOperator, 2);
nPort = numel(coilCurrent);
weighted = reshape(reshape(forceOperator, 3 * nState, nPort) * coilCurrent(:), 3, nState);
force_N = weighted * modalState(:);
end
