function result = canonicalBodyHamiltonianJet(coefficients, magneticRigidityTM, options)
%CANONICALBODYHAMILTONIANJET Build the native fourth-degree body Hamiltonian.
%   RESULT = radia.beam.canonicalBodyHamiltonianJet(COEFFICIENTS, BRHO)
%   accepts dipole, normal/skew quadrupole, normal/skew sextupole, and
%   normal/skew octupole coefficients. It returns symmetric H2/H3/H4 and
%   the canonical dynamics A/F2/F3 in coordinates
%   (x,px/p0,y,py/p0,ell,delta). The longitudinal Poisson sign is -1.
arguments
    coefficients {mustBeNumeric,mustBeReal,mustBeFinite}
    magneticRigidityTM (1,1) double {mustBeFinite,mustBeNonzero}
    options.CurvatureSign (1,1) double {mustBeFinite} = 1
    options.GradientSign (1,1) double {mustBeFinite} = 1
    options.ReferenceBeta (1,1) double {mustBeFinite,mustBePositive} = 1
end
values = double(coefficients(:).');
if numel(values) ~= 7
    error("radia:beam:InvalidShape", ...
        "coefficients must contain exactly seven entries.");
end
if options.ReferenceBeta > 1
    error("radia:beam:InvalidReferenceBeta", ...
        "ReferenceBeta must be in (0,1].");
end
config = struct( ...
    schema='radia.beam.canonical-hamiltonian-jet.v1', ...
    coefficients=values, ...
    magnetic_rigidity_t_m=magneticRigidityTM, ...
    curvature_sign=options.CurvatureSign, ...
    gradient_sign=options.GradientSign, ...
    reference_beta=options.ReferenceBeta);
result = radia.internal.callMex( ...
    "beam.hamiltonian.canonical_body_jet",config);
end
