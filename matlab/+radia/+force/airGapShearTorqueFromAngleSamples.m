function summary = airGapShearTorqueFromAngleSamples( ...
    angles_rad, magneticFluxDensityRadial_T, magneticFluxDensityTangential_T, ...
    radius_m, axialLength_m, periodic, period_rad, permeability_H_per_m)
%AIRGAPSHEARTORQUEFROMANGLESAMPLES Integrate sampled cylindrical-gap torque.

if nargin < 5 || isempty(axialLength_m), axialLength_m = 1.0; end
if nargin < 6 || isempty(periodic), periodic = true; end
if nargin < 7 || isempty(period_rad), period_rad = 2*pi; end
if nargin < 8 || isempty(permeability_H_per_m), permeability_H_per_m = 4*pi*1e-7; end
angles = localTable(angles_rad, "angles_rad");
radial = localTable(magneticFluxDensityRadial_T, "magneticFluxDensityRadial_T");
tangential = localTable(magneticFluxDensityTangential_T, "magneticFluxDensityTangential_T");
if numel(angles) ~= numel(radial) || numel(angles) ~= numel(tangential)
    error("radia:force:SampleCount", "angles and air-gap fields must have the same length");
end
if any(diff(angles) <= 0)
    error("radia:force:Angles", "angles_rad must be strictly increasing");
end
localGeometry(radius_m, axialLength_m, permeability_H_per_m);
shear = radial .* tangential / permeability_H_per_m;
[integral, integratedAngle, rows] = localIntegrate(angles, shear, logical(periodic), period_rad);
force = radius_m * axialLength_m * integral;
torque = radius_m * force;
if axialLength_m > 0
    torquePerLength = torque/axialLength_m;
else
    torquePerLength = inf;
end
for index = 1:numel(rows)
    rows(index).tangential_force_N = rows(index).shear_average_Pa * ...
        radius_m * axialLength_m * rows(index).angle_width_rad;
    rows(index).torque_Nm = rows(index).tangential_force_N * radius_m;
end
summary = struct( ...
    "n_samples", numel(angles), "n_segments", numel(rows), ...
    "periodic", logical(periodic), "period_rad", period_rad, ...
    "radius_m", radius_m, "axial_length_m", axialLength_m, ...
    "permeability_H_per_m", permeability_H_per_m, ...
    "integrated_angle_rad", integratedAngle, ...
    "integral_shear_dtheta_Pa_rad", integral, ...
    "average_shear_stress_Pa", integral/integratedAngle, ...
    "tangential_force_N", force, "torque_Nm", torque, ...
    "torque_per_axial_length_N", torquePerLength, "rows", rows);
end

function values = localTable(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || numel(value) < 2 || any(~isfinite(value), "all")
    error("radia:force:Table", "%s must be a finite real vector with at least two samples", name);
end
values = double(value(:));
end

function localGeometry(radius, lengthValue, permeability)
if any(~isfinite([radius, lengthValue])) || any([radius, lengthValue] < 0)
    error("radia:force:Geometry", "radius_m and axialLength_m must be finite and nonnegative");
end
if ~isscalar(permeability) || ~isfinite(permeability) || permeability <= 0
    error("radia:force:Permeability", "permeability_H_per_m must be finite and positive");
end
end

function [integral, integratedAngle, rows] = localIntegrate(angles, shear, periodic, period)
if ~isscalar(period) || ~isfinite(period) || period <= 0
    error("radia:force:Period", "period_rad must be finite and positive");
end
if periodic && angles(end) - angles(1) >= period
    error("radia:force:Period", "periodic angles must omit the duplicate endpoint");
end
count = numel(angles);
segmentCount = count - 1 + double(periodic);
rows = repmat(struct("segment_index", 0, "angle_start_rad", 0, ...
    "angle_end_rad", 0, "angle_width_rad", 0, "shear_start_Pa", 0, ...
    "shear_end_Pa", 0, "shear_average_Pa", 0, ...
    "tangential_force_N", 0, "torque_Nm", 0), segmentCount, 1);
integral = 0;
integratedAngle = 0;
for index = 1:segmentCount
    nextIndex = mod(index, count) + 1;
    angleEnd = angles(nextIndex);
    if periodic && nextIndex == 1, angleEnd = angleEnd + period; end
    width = angleEnd - angles(index);
    average = 0.5 * (shear(index) + shear(nextIndex));
    integral = integral + average * width;
    integratedAngle = integratedAngle + width;
    rows(index).segment_index = index;
    rows(index).angle_start_rad = angles(index);
    rows(index).angle_end_rad = angleEnd;
    rows(index).angle_width_rad = width;
    rows(index).shear_start_Pa = shear(index);
    rows(index).shear_end_Pa = shear(nextIndex);
    rows(index).shear_average_Pa = average;
end
end
