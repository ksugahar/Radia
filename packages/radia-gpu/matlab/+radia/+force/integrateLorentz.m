function force_N = integrateLorentz(currentDensity_A_per_m2, ...
    magneticFluxDensity_T, volumeWeights_m3)
%INTEGRATELORENTZ Integrate J cross B over physical volume weights.
%   The return value is the Cartesian 1-by-3 force vector in newtons.

forceDensity = radia.force.lorentzDensity( ...
    currentDensity_A_per_m2, magneticFluxDensity_T);
weights = localWeights(volumeWeights_m3, size(forceDensity, 1), ...
    "volumeWeights_m3");
force_N = sum(forceDensity .* weights, 1);
end

function weights = localWeights(value, count, name)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all")
    error("radia:force:Weights", "%s must contain finite real values", name);
end
weights = double(value(:));
if any(weights < 0)
    error("radia:force:Weights", "%s must be nonnegative", name);
end
if isscalar(weights)
    weights = repmat(weights, count, 1);
elseif numel(weights) ~= count
    error("radia:force:WeightCount", ...
        "%s must be scalar or contain one value per sample", name);
end
end
