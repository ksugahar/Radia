function force_N = evaluateHCurlEddyFosterForce(model, modalCoefficients, coilCurrent)
%EVALUATEHCURLEDDYFOSTERFORCE Evaluate the modal reduced Lorentz force.
%   K(k,j,b) is the reduced force operator contracted onto the Foster modes
%   (component k, mode j, coil port b).  The phasor time average is
%       F_k = 0.5*real(sum_j sum_b K(k,j,b)*z_j*conj(i_b)).
%   These are phasors.  The Simulink block evaluates the instantaneous
%   product K*z(t)*i(t) of real signals, whose cycle average is this value.

arguments
    model (1,1) struct
    modalCoefficients double {mustBeFinite}
    coilCurrent double {mustBeFinite}
end

nState = model.state_order;
nPort = model.port_count;
K = double(model.force_operator);
if size(K, 1) ~= 3 || size(K, 2) ~= nState || numel(K) ~= 3 * nState * nPort
    error("radia:simulink:HCurlFosterForce", ...
        "force_operator must have size [3, n_state, port_count].");
end
z = double(modalCoefficients(:));
if numel(z) ~= nState
    error("radia:simulink:HCurlFosterForce", ...
        "modalCoefficients must contain n_state values.");
end
i = double(coilCurrent(:));
if isscalar(i)
    i = repmat(i, nPort, 1);
elseif numel(i) ~= nPort
    error("radia:simulink:HCurlFosterForce", ...
        "coilCurrent must be scalar or contain port_count values.");
end
weighted = reshape(reshape(K, 3 * nState, nPort) * conj(i), 3, nState);
force_N = 0.5 * real(weighted * z);
end
