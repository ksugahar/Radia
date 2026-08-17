function result = canonicalBodyHamiltonianJet(coefficients, magneticRigidityTM, options)
%CANONICALBODYHAMILTONIANJET Build the native fifth-degree body Hamiltonian.
%   RESULT = radia.beam.canonicalBodyHamiltonianJet(COEFFICIENTS, BRHO)
%   accepts dipole, normal/skew quadrupole, normal/skew sextupole, and
%   normal/skew octupole and optional normal/skew decapole coefficients. It
%   returns symmetric H2/H3/H4/H5 and canonical dynamics A/F2/F3/F4 in coordinates
%   (x,px/p0,y,py/p0,ell,delta). The longitudinal Poisson sign is -1.
arguments
    coefficients {mustBeNumeric,mustBeReal,mustBeFinite}
    magneticRigidityTM (1,1) double {mustBeFinite,mustBeNonzero}
    options.CurvatureSign (1,1) double {mustBeFinite} = 1
    options.GradientSign (1,1) double {mustBeFinite} = 1
    options.ReferenceBeta (1,1) double {mustBeFinite,mustBePositive} = 1
    options.ReferenceCurvaturePerM {mustBeNumeric,mustBeReal} = []
end
values = double(coefficients(:).');
if ~ismember(numel(values),[7,9])
    error("radia:beam:InvalidShape", ...
        "coefficients must contain seven or nine entries.");
end
if options.ReferenceBeta > 1
    error("radia:beam:InvalidReferenceBeta", ...
        "ReferenceBeta must be in (0,1].");
end
if ~(isempty(options.ReferenceCurvaturePerM) || ...
        (isscalar(options.ReferenceCurvaturePerM) && ...
        isfinite(options.ReferenceCurvaturePerM)))
    error("radia:beam:InvalidReferenceCurvature", ...
        "ReferenceCurvaturePerM must be empty or one finite scalar.");
end
config = struct( ...
    schema='radia.beam.canonical-hamiltonian-jet.v1', ...
    coefficients=values, ...
    magnetic_rigidity_t_m=magneticRigidityTM, ...
    curvature_sign=options.CurvatureSign, ...
    gradient_sign=options.GradientSign, ...
    reference_beta=options.ReferenceBeta);
if ~isempty(options.ReferenceCurvaturePerM)
    config.reference_curvature_per_m = ...
        double(options.ReferenceCurvaturePerM);
end
result = radia.internal.callMex( ...
    "beam.hamiltonian.canonical_body_jet",config);
end
