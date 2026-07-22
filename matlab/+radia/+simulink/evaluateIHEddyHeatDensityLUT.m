function heatDensity_W_per_m3 = evaluateIHEddyHeatDensityLUT( ...
        lut, coilCurrentRms_A, rotationAngle_rad)
%EVALUATEIHEDDYHEATDENSITYLUT Evaluate a periodic IH loss-density map.
%   Q = radia.simulink.evaluateIHEddyHeatDensityLUT(LUT,I,THETA) returns
%   one row per input sample and one column per heated region. Scalar input
%   expands against the other input. Current sign is ignored because LUT
%   current is the RMS carrier-current amplitude; angles wrap periodically.

arguments
    lut (1,1) struct
    coilCurrentRms_A double {mustBeFinite}
    rotationAngle_rad double {mustBeFinite}
end

if ~isfield(lut, "schema") || ...
        lut.schema ~= "radia.ih.simulink.eddy_heat_density_lut.v1"
    error("radia:simulink:InvalidEddyLUT", ...
        "unsupported IH eddy heat-density LUT schema.");
end

current = coilCurrentRms_A(:);
theta = rotationAngle_rad(:);
if isscalar(current) && ~isscalar(theta)
    current = repmat(current, size(theta));
elseif isscalar(theta) && ~isscalar(current)
    theta = repmat(theta, size(current));
elseif numel(current) ~= numel(theta)
    error("radia:simulink:EddyLUTInputSize", ...
        "current and rotation angle must be scalar or have the same number of samples.");
end

theta = mod(theta - lut.angle_origin_rad, lut.angle_period_rad) + ...
    lut.angle_origin_rad;
current = abs(current);
currentGrid = lut.coil_current_rms_breakpoints_A;
current = min(max(current, currentGrid(1)), currentGrid(end));

nSample = numel(current);
nRegion = lut.region_count;
heatDensity_W_per_m3 = zeros(nSample, nRegion);
for region = 1:nRegion
    interpolant = griddedInterpolant( ...
        {lut.rotation_angle_breakpoints_rad, currentGrid}, ...
        lut.table_heat_density_W_per_m3(:, :, region), ...
        "linear", "nearest");
    heatDensity_W_per_m3(:, region) = interpolant(theta, current);
end
end
