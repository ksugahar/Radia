function model = hcurl_eddy_foster_model(vol_path, order, ports, steps, options)
%HCURL_EDDY_FOSTER_MODEL Build a native HCurl diffusion Foster model from a VOL mesh.
%   MODEL = radia.ngsolve.hcurl_eddy_foster_model(VOL_PATH, ORDER, PORTS, STEPS)
%   builds the high-order HCurl response basis in the C++ MEX gateway and
%   projects the NGSolve mass and curl-curl operators into that basis:
%       M_r = V' * M * V,  K_r = V' * K * V,  P_r = V' * PORTS,
%   followed by the local HCurl diffusion convention
%       R = Reluctivity * K_r,  L = Conductivity * M_r.
%   The reduced pair is diagonalised by R*w = lambda*L*w with W'*L*W = I, and
%   the result is a Foster state-space contract (makeHCurlEddyFosterModel)
%   with decay rates lambda and modal input W'*P_r.  This is a Python-free
%   local FE path; it does not replace the separate VIM/BEM Laplace
%   inductance or a frequency-dependent SIBC rationalization.

arguments
    vol_path (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    ports double {mustBeReal, mustBeFinite, mustBeNonempty}
    steps (1,1) double {mustBeInteger, mustBePositive}
    options.NoGrads (1,1) logical = true
    options.Rtol (1,1) double {mustBePositive, mustBeFinite} = 1.0e-12
    options.Conductivity (1,1) double {mustBePositive, mustBeFinite} = 1.0
    options.Reluctivity (1,1) double {mustBePositive, mustBeFinite} = 1.0
    options.SampleTime_s (1,1) double {mustBePositive, mustBeFinite} = 1.0e-5
    options.PassivityTolerance (1,1) double {mustBeNonnegative} = 1.0e-10
    options.InitialState double = []
end

if ~ismatrix(ports)
    error("radia:ngsolve:HCurlPorts", ...
        "ports must be a two-dimensional ndof-by-nports real matrix.");
end

basis = radia.ngsolve.hcurl_eddy_native_basis( ...
    vol_path, order, ports, steps, ...
    no_grads=options.NoGrads, rtol=options.Rtol);
if basis.rank < 1
    error("radia:ngsolve:HCurlReduction", ...
        "the native HCurl response reduction returned rank zero.");
end

resistance = options.Reluctivity * basis.curlcurl_gram;
inductance = options.Conductivity * basis.mass_gram;
resistance = 0.5 * (resistance + resistance.');
inductance = 0.5 * (inductance + inductance.');
[~, notPositive] = chol(inductance);
if notPositive
    error("radia:ngsolve:HCurlInductance", ...
        "the reduced inductance Conductivity*V'*M*V is not positive definite.");
end
[modes, rates] = eig(resistance, inductance, "chol", "vector");
modes = modes ./ sqrt(sum(modes .* (inductance * modes), 1));
[rates, order_] = sort(rates);
modes = modes(:, order_);

model = radia.simulink.makeHCurlEddyFosterModel(rates, modes.' * basis.port_rhs, ...
    SampleTime_s=options.SampleTime_s, ...
    PassivityTolerance=options.PassivityTolerance, ...
    InitialState=options.InitialState);
model.assembly_schema = "radia.hcurl.eddy_foster.native_diffusion.v1";
model.parent_space = "NGSolve HCurl";
model.vol_path = vol_path;
model.parent_order = order;
model.krylov_steps = steps;
model.conductivity = options.Conductivity;
model.reluctivity = options.Reluctivity;
model.resistance = resistance;
model.inductance = inductance;
model.modes = modes;
model.native_basis = basis;
model.projection = string(basis.projection);
end
