function force_N = integrateTimeAverageMaxwellSurface( ...
    magneticFluxDensityPhasor_T, outwardNormal, areaWeights_m2, ...
    permeability_H_per_m, amplitude)
%INTEGRATETIMEAVERAGEMAXWELLSURFACE Integrate cycle-averaged air traction.

if nargin < 4
    permeability_H_per_m = [];
end
if nargin < 5
    amplitude = [];
end
traction = radia.force.timeAverageMaxwellTractionAir( ...
    magneticFluxDensityPhasor_T, outwardNormal, permeability_H_per_m, amplitude);
weights = localWeights(areaWeights_m2, size(traction, 1));
force_N = sum(traction .* weights, 1);
end

function weights = localWeights(value, count)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all")
    error("radia:force:Weights", "areaWeights_m2 must contain finite real values");
end
weights = double(value(:));
if any(weights < 0)
    error("radia:force:Weights", "areaWeights_m2 must be nonnegative");
end
if isscalar(weights)
    weights = repmat(weights, count, 1);
elseif numel(weights) ~= count
    error("radia:force:WeightCount", "areaWeights_m2 must be scalar or contain one value per sample");
end
end
