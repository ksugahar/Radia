function model = hcurl_eddy_cln_model(vol_path, order, ports, steps, options)
%HCURL_EDDY_CLN_MODEL Build a native HCurl diffusion CLN model from a VOL mesh.
%   MODEL = radia.ngsolve.hcurl_eddy_cln_model(VOL_PATH, ORDER, PORTS, STEPS)
%   builds the high-order HCurl response basis in the C++ MEX gateway and
%   projects the NGSolve mass and curl-curl operators into that basis.
%
%   The native projection is
%       M_r = V' * M * V,  K_r = V' * K * V,  P_r = V' * PORTS,
%   followed by the local HCurl diffusion convention
%       R = Reluctivity * K_r,  L = Conductivity * M_r.
%   This is a Python-free local FE CLN path.  It does not claim to replace
%   the separate VIM/BEM Laplace inductance or a frequency-dependent SIBC
%   rationalization; those remain explicit model stages.

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

if ndims(ports) ~= 2
    error("radia:ngsolve:HCurlPorts", ...
        "ports must be a two-dimensional ndof-by-nports real matrix.");
end

basis = radia.ngsolve.hcurl_eddy_cln_native_basis( ...
    vol_path, order, ports, steps, ...
    no_grads=options.NoGrads, rtol=options.Rtol);
if basis.rank < 1
    error("radia:ngsolve:HCurlReduction", ...
        "the native HCurl response reduction returned rank zero.");
end

resistance = options.Reluctivity * basis.curlcurl_gram;
inductance = options.Conductivity * basis.mass_gram;
portRHS = basis.port_rhs;
model = radia.simulink.makeHCurlEddyCLNModel( ...
    resistance, inductance, portRHS, ...
    SampleTime_s=options.SampleTime_s, ...
    PassivityTolerance=options.PassivityTolerance, ...
    InitialState=options.InitialState);

model.assembly_schema = "radia.hcurl.eddy_cln.native_diffusion.v1";
model.parent_space = "NGSolve HCurl";
model.vol_path = vol_path;
model.parent_order = order;
model.krylov_steps = steps;
model.conductivity = options.Conductivity;
model.reluctivity = options.Reluctivity;
model.native_basis = basis;
model.projection = string(basis.projection);
end
