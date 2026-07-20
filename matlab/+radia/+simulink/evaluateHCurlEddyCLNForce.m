function force_N = evaluateHCurlEddyCLNForce(model, coefficients, coilCurrent)
%EVALUATEHCURLEDDYCLNFORCE Evaluate the reduced Lorentz force operator.
%   The force operator is exported by NGSolve as K(k,a,b), where k is the
%   Cartesian component, a is the reduced current mode, and b is the coil
%   port.  The phasor time-average is
%       F_k = 0.5*real(sum(K(k,a,b)*c_a*conj(i_b))).

arguments
    model (1,1) struct
    coefficients double {mustBeFinite}
    coilCurrent double {mustBeFinite}
end

if ~isfield(model, "force_operator")
    error("radia:simulink:HCurlCLNForce", ...
        "model does not contain a force_operator.");
end
K = double(model.force_operator);
nState = model.state_order;
nPort = model.port_count;
if numel(K) ~= 3 * nState * nPort || ...
        size(K, 1) ~= 3 || size(K, 2) ~= nState
    error("radia:simulink:HCurlCLNForce", ...
        "force_operator must have size [3, n_state, port_count].");
end
K = reshape(K, [3, nState, nPort]);
c = double(coefficients(:));
if numel(c) ~= nState
    error("radia:simulink:HCurlCLNForce", ...
        "coefficients must contain n_state values.");
end
i = double(coilCurrent(:));
if numel(i) == 1 && nPort > 1
    i = repmat(i, nPort, 1);
end
if numel(i) ~= nPort
    error("radia:simulink:HCurlCLNForce", ...
        "coilCurrent must contain port_count values.");
end
force_N = zeros(3, 1);
for k = 1:3
    for a = 1:nState
        for b = 1:nPort
            force_N(k) = force_N(k) + K(k, a, b) * c(a) * conj(i(b));
        end
    end
end
force_N = 0.5 * real(force_N);
end
