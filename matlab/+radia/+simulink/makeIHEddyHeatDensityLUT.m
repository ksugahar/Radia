function lut = makeIHEddyHeatDensityLUT(rotationAngle_rad, ...
        coilCurrentRms_A, heatDensity_W_per_m3, options)
%MAKEIHEDDYHEATDENSITYLUT Create a periodic IH eddy-loss density map.
%   LUT = radia.simulink.makeIHEddyHeatDensityLUT(THETA,I,Q) validates a
%   Radia-generated table Q(theta,current,region). THETA contains samples
%   from one mechanical period without requiring the repeated endpoint.
%   I contains nonnegative RMS carrier-current amplitudes and must start at
%   zero. Q may be angle-by-current for one heated region or
%   angle-by-current-by-region for a partitioned workpiece.
%
%   The returned table includes the periodic endpoint and is suitable for
%   native Simulink 2-D Lookup Table blocks. RegionVolumes_m3 defines the
%   density-to-total-power reduction used by the lumped Thermal subsystem.

arguments
    rotationAngle_rad (:,1) double {mustBeFinite}
    coilCurrentRms_A (:,1) double {mustBeFinite}
    heatDensity_W_per_m3 double
    options.AnglePeriod_rad (1,1) double {mustBePositive} = 2 * pi
    options.RegionVolumes_m3 (:,1) double {mustBePositive} = 1.0
    options.RegionNames (:,1) string = string.empty(0, 1)
    options.CarrierFrequency_Hz (1,1) double {mustBeNonnegative} = 0.0
    options.Source (1,1) string = "Radia VIM/FEM/SIBC/ESIM"
end

theta = rotationAngle_rad(:);
current = coilCurrentRms_A(:);
nAngle = numel(theta);
nCurrent = numel(current);
nRegion = numel(options.RegionVolumes_m3);
if nAngle < 2 || any(diff(theta) <= 0)
    error("radia:simulink:EddyAngleGrid", ...
        "rotationAngle_rad must contain at least two increasing samples.");
end
if nCurrent < 2 || any(current < 0) || any(diff(current) <= 0)
    error("radia:simulink:EddyCurrentGrid", ...
        "coilCurrentRms_A must contain at least two increasing nonnegative samples.");
end
currentTolerance = 100 * eps(max(1, current(end)));
if abs(current(1)) > currentTolerance
    error("radia:simulink:EddyCurrentOrigin", ...
        "coilCurrentRms_A must start at zero so zero-current heating is explicit.");
end

actualSize = size(heatDensity_W_per_m3);
actualSize(end + 1:3) = 1;
if actualSize(1) ~= nAngle || actualSize(2) ~= nCurrent || ...
        actualSize(3) ~= nRegion
    error("radia:simulink:EddyHeatDensityShape", ...
        "heatDensity_W_per_m3 must have size [%d %d %d].", ...
        nAngle, nCurrent, nRegion);
end
table = reshape(heatDensity_W_per_m3, nAngle, nCurrent, nRegion);
if any(~isfinite(table(:))) || any(table(:) < 0)
    error("radia:simulink:EddyHeatDensity", ...
        "heatDensity_W_per_m3 must be finite and nonnegative.");
end
zeroTolerance = 1.0e3 * eps(max(1, max(table(:))));
if any(abs(table(:, 1, :)) > zeroTolerance, "all")
    error("radia:simulink:EddyZeroCurrent", ...
        "the zero-current heat-density column must be zero.");
end

period = options.AnglePeriod_rad;
span = theta(end) - theta(1);
periodTolerance = 1.0e3 * eps(max(1, period));
if span > period + periodTolerance
    error("radia:simulink:EddyAnglePeriod", ...
        "rotation-angle samples must not span more than AnglePeriod_rad.");
end
if abs(span - period) <= periodTolerance
    endpointTolerance = 1.0e3 * eps(max(1, max(table(:))));
    if any(abs(table(end, :, :) - table(1, :, :)) > endpointTolerance, "all")
        error("radia:simulink:EddyPeriodicEndpoint", ...
            "heat density at the repeated period endpoint must equal the first row.");
    end
else
    theta(end + 1, 1) = theta(1) + period;
    table(end + 1, :, :) = table(1, :, :);
end

if options.RegionNames.isempty()
    regionNames = "region" + string((1:nRegion).');
elseif numel(options.RegionNames) ~= nRegion
    error("radia:simulink:EddyRegionNames", ...
        "RegionNames must match RegionVolumes_m3.");
else
    regionNames = options.RegionNames(:);
end

lut = struct( ...
    "schema", "radia.ih.simulink.eddy_heat_density_lut.v1", ...
    "rotation_angle_breakpoints_rad", theta, ...
    "coil_current_rms_breakpoints_A", current, ...
    "table_heat_density_W_per_m3", table, ...
    "angle_period_rad", period, ...
    "angle_origin_rad", theta(1), ...
    "region_names", regionNames, ...
    "region_volumes_m3", options.RegionVolumes_m3(:), ...
    "region_count", nRegion, ...
    "carrier_frequency_Hz", options.CarrierFrequency_Hz, ...
    "input_names", ["coil_current_rms_A"; "rotation_angle_rad"], ...
    "output_name", "heat_density_W_per_m3", ...
    "source", options.Source, ...
    "current_extrapolation", "clip", ...
    "angle_extrapolation", "periodic");
end
