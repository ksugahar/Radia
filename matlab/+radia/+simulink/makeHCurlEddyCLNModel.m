function model = makeHCurlEddyCLNModel(resistance, inductance, portRHS, options)
%MAKEHCURLEDDYCLNMODEL Create a MATLAB contract for a reduced HCurl-VIM.
%   The reduced model is
%       (R + s*L)c = -s*P*i,
%   where P is the vector-potential port matrix.  The time-domain input
%   convention is u = -di/dt, giving c_dot = A*c + B*u and y = P'*c.
%   This is the numeric MATLAB boundary for a p=6 HCurl Eddy Bubble/CLN
%   model assembled by NGSolve or another trusted frontend.

arguments
    resistance double
    inductance double
    portRHS double
    options.SampleTime_s (1,1) double {mustBePositive} = 1.0e-5
    options.PassivityTolerance (1,1) double {mustBeNonnegative} = 1.0e-10
    options.InitialState double = []
end

R = double(resistance);
L = double(inductance);
P = double(portRHS);
if ndims(R) ~= 2 || size(R, 1) ~= size(R, 2) || isempty(R)
    error("radia:simulink:HCurlCLNMatrix", ...
        "resistance must be a non-empty square matrix.");
end
n = size(R, 1);
if ~isequal(size(L), [n, n])
    error("radia:simulink:HCurlCLNMatrix", ...
        "inductance must match resistance dimensions.");
end
if ndims(P) ~= 2 || size(P, 1) ~= n || size(P, 2) < 1
    error("radia:simulink:HCurlCLNMatrix", ...
        "portRHS must have n_state rows and at least one port column.");
end
if any(~isfinite(R), "all") || any(~isfinite(L), "all") || any(~isfinite(P), "all")
    error("radia:simulink:HCurlCLNMatrix", "HCurl CLN matrices must be finite.");
end

R = 0.5 * (R + R.');
L = 0.5 * (L + L.');
if min(eig(R)) < -options.PassivityTolerance || ...
        min(eig(L)) < -options.PassivityTolerance
    error("radia:simulink:HCurlCLNPassive", ...
        "resistance and inductance must be positive semidefinite.");
end
if rcond(L) <= eps
    error("radia:simulink:HCurlCLNMatrix", ...
        "inductance is singular to working precision.");
end

A = -(L \ R);
B = L \ P;
C = P.';
D = zeros(size(C, 1), size(B, 2));
if isempty(options.InitialState)
    x0 = zeros(n, 1);
else
    x0 = double(options.InitialState(:));
    if numel(x0) ~= n || any(~isfinite(x0))
        error("radia:simulink:HCurlCLNState", ...
            "InitialState must contain n_state finite values.");
    end
end

% Exact zero-order-hold discretization keeps the Simulink block faithful to
% the continuous reduced model without requiring Control System Toolbox.
augmented = [A, B; zeros(size(B, 2), n + size(B, 2))];
discrete = expm(augmented * options.SampleTime_s);
Ad = discrete(1:n, 1:n);
Bd = discrete(1:n, n + 1:end);

model = struct( ...
    "schema", "radia.hcurl.eddy_cln.state_space.v1", ...
    "resistance", R, "inductance", L, "port_rhs", P, ...
    "A", A, "B", B, "C", C, "D", D, ...
    "Ad", Ad, "Bd", Bd, "Cd", C, "Dd", D, ...
    "x0", x0, "sample_time_s", options.SampleTime_s, ...
    "state_order", n, "port_count", size(P, 2), ...
    "passive", true, ...
    "input_convention", "u=-d(coil_current)/dt", ...
    "has_sibc_termination", false, ...
    "notes", "Reduced HCurl Eddy Bubble/CLN model; SIBC must be rationalized before state-space export.");
end
