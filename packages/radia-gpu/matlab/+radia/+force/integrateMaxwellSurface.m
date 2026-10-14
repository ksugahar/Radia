function force_N = integrateMaxwellSurface(magneticFluxDensity_T, ...
    outwardNormal, areaWeights_m2, permeability_H_per_m)
%INTEGRATEMAXWELLSURFACE Integrate air-side Maxwell traction over a surface.
%   The surface must be closed around the requested body and remain in air.
%   AREAWEIGHTS_M2 includes the physical surface Jacobian and quadrature weight.

arguments
    magneticFluxDensity_T {mustBeNumeric,mustBeReal,mustBeFinite}
    outwardNormal {mustBeNumeric,mustBeReal,mustBeFinite}
    areaWeights_m2 {mustBeNumeric,mustBeReal,mustBeFinite}
    permeability_H_per_m (1,1) double {mustBeFinite,mustBePositive} = 4*pi*1e-7
end

traction = radia.force.maxwellTractionAir( ...
    magneticFluxDensity_T, outwardNormal, permeability_H_per_m);
weights = double(areaWeights_m2(:));
if any(weights < 0)
    error("radia:force:Weights", "areaWeights_m2 must be nonnegative");
end
if isscalar(weights)
    weights = repmat(weights, size(traction, 1), 1);
elseif numel(weights) ~= size(traction, 1)
    error("radia:force:WeightCount", ...
        "areaWeights_m2 must be scalar or contain one value per sample");
end
force_N = sum(traction .* weights, 1);
end
