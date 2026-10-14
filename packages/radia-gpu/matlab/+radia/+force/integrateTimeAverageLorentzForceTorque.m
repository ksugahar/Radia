function [force_N, torque_Nm] = integrateTimeAverageLorentzForceTorque( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, ...
    volumeWeights_m3, samplePoints_m, pivot_m, amplitude)
%INTEGRATETIMEAVERAGELORENTZFORCETORQUE Integrate phasor resultants.

if nargin < 5 || isempty(pivot_m)
    pivot_m = [0, 0, 0];
end
if nargin < 6 || isempty(amplitude)
    amplitude = "peak";
end
density = radia.force.timeAverageLorentzDensity( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, amplitude);
[points, pivot, weights] = localGeometry( ...
    samplePoints_m, pivot_m, volumeWeights_m3, size(density, 1));
force_N = sum(density .* weights, 1);
torque_Nm = sum(cross(points - pivot, density, 2) .* weights, 1);
end

function [points, pivot, weights] = localGeometry(pointsValue, pivotValue, weightsValue, count)
if ~isnumeric(pointsValue) || ~isreal(pointsValue) || any(~isfinite(pointsValue), "all")
    error("radia:force:SamplePoints", "samplePoints_m must contain finite real values");
end
points = double(pointsValue);
if isvector(points) && numel(points) == 3
    points = reshape(points, 1, 3);
end
if ndims(points) ~= 2 || size(points, 2) ~= 3 || size(points, 1) ~= count %#ok<ISMAT>
    error("radia:force:SamplePointCount", "samplePoints_m must be N-by-3 with one row per sample");
end
if ~isnumeric(pivotValue) || ~isreal(pivotValue) || numel(pivotValue) ~= 3 || any(~isfinite(pivotValue), "all")
    error("radia:force:Pivot", "pivot_m must be one finite real three-vector");
end
pivot = reshape(double(pivotValue), 1, 3);
weights = double(weightsValue(:));
if ~isnumeric(weightsValue) || ~isreal(weightsValue) || any(~isfinite(weights), "all") || any(weights < 0)
    error("radia:force:Weights", "volumeWeights_m3 must contain finite nonnegative real values");
end
if isscalar(weights)
    weights = repmat(weights, count, 1);
elseif numel(weights) ~= count
    error("radia:force:WeightCount", "volumeWeights_m3 must be scalar or contain one value per sample");
end
end
