function force_N = integrateTimeAverageLorentz( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, ...
    volumeWeights_m3, amplitude)
%INTEGRATETIMEAVERAGELORENTZ Integrate cycle-averaged Lorentz force.

if nargin < 4
    amplitude = "peak";
end
density = radia.force.timeAverageLorentzDensity( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, amplitude);
weights = localWeights(volumeWeights_m3, size(density, 1));
force_N = sum(density .* weights, 1);
end

function weights = localWeights(value, count)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all")
    error("radia:force:Weights", "volumeWeights_m3 must contain finite real values");
end
weights = double(value(:));
if any(weights < 0)
    error("radia:force:Weights", "volumeWeights_m3 must be nonnegative");
end
if isscalar(weights)
    weights = repmat(weights, count, 1);
elseif numel(weights) ~= count
    error("radia:force:WeightCount", "volumeWeights_m3 must be scalar or contain one value per sample");
end
end
