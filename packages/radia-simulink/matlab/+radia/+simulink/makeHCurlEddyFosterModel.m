function model = makeHCurlEddyFosterModel(decayRates, modalPortRHS, options)
%MAKEHCURLEDDYFOSTERMODEL Create the MATLAB contract for a Foster HCurl model.
%   The reduced HCurl eddy-current model (R + s*L)c = -s*P*i is diagonalised
%   by R*v = lambda*L*v with V'*L*V = I.  In modal coordinates z (c = V*z)
%   and for the input u = -di/dt,
%       z_dot = -diag(lambda)*z + Bm*u,   y = Bm.'*z,   Bm = V.'*P.
%   radia.vim.ExportHCurlEddyFosterJSON writes lambda and Bm.  The optional
%   force operator is the reduced Lorentz operator contracted onto the
%   modes, K(k,j,b).

arguments
    decayRates double
    modalPortRHS double
    options.SampleTime_s (1,1) double {mustBeFinite, mustBePositive} = 1.0e-5
    options.ForceOperator double = []
    options.InitialState double = []
    options.PassivityTolerance (1,1) double {mustBeNonnegative} = 1.0e-10
end

lambda = double(decayRates(:));
Bm = double(modalPortRHS);
n = numel(lambda);
if n < 1 || any(~isfinite(lambda))
    error("radia:simulink:HCurlFosterModel", ...
        "decayRates must be a non-empty finite vector.");
end
scale = max([1.0; abs(lambda)]);
if any(lambda < -options.PassivityTolerance * scale)
    error("radia:simulink:HCurlFosterPassive", ...
        "decay rates must be nonnegative for a passive Foster model.");
end
lambda = max(lambda, 0.0);
if ~ismatrix(Bm) || size(Bm, 1) ~= n || size(Bm, 2) < 1 || ...
        any(~isfinite(Bm), "all")
    error("radia:simulink:HCurlFosterModel", ...
        "modalPortRHS must be a finite n_state-by-port_count matrix.");
end
nPort = size(Bm, 2);

if isempty(options.ForceOperator)
    K = zeros(3, n, nPort);
else
    K = double(options.ForceOperator);
    if numel(K) ~= 3 * n * nPort || size(K, 1) ~= 3 || size(K, 2) ~= n || ...
            any(~isfinite(K), "all")
        error("radia:simulink:HCurlFosterForce", ...
            "ForceOperator must be a finite [3, n_state, port_count] array.");
    end
    K = reshape(K, [3, n, nPort]);
end
if isempty(options.InitialState)
    x0 = zeros(n, 1);
else
    x0 = double(options.InitialState(:));
    if numel(x0) ~= n || any(~isfinite(x0))
        error("radia:simulink:HCurlFosterState", ...
            "InitialState must contain n_state finite values.");
    end
end

[adDiagonal, bdGain] = radia.internal.hcurlEddyFosterZOH(lambda, options.SampleTime_s);
C = Bm.';
D = zeros(nPort, nPort);
model = struct( ...
    "schema", "radia.hcurl.eddy_foster.state_space.v1", ...
    "decay_rates", lambda, "modal_port_rhs", Bm, ...
    "A", -diag(lambda), "B", Bm, "C", C, "D", D, ...
    "Ad", diag(adDiagonal), "Bd", bdGain .* Bm, "Cd", C, "Dd", D, ...
    "force_operator", K, ...
    "x0", x0, "sample_time_s", options.SampleTime_s, ...
    "state_order", n, "port_count", nPort, ...
    "passive", true, ...
    "input_convention", "u=-d(coil_current)/dt", ...
    "state_convention", "z_dot=-diag(decay_rates)*z+modal_port_rhs*u, y=modal_port_rhs.'*z", ...
    "has_sibc_termination", false);
end
