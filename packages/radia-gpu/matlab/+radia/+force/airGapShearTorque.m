function torque_Nm = airGapShearTorque( ...
    magneticFluxDensityRadial_T, magneticFluxDensityTangential_T, ...
    radius_m, axialLength_m, angle_rad, permeability_H_per_m)
%AIRGAPSHEARTORQUE Uniform cylindrical-gap Maxwell shear torque in N m.

if nargin < 4 || isempty(axialLength_m)
    axialLength_m = 1.0;
end
if nargin < 5 || isempty(angle_rad)
    angle_rad = 2*pi;
end
if nargin < 6
    permeability_H_per_m = [];
end
values = [radius_m, axialLength_m, angle_rad];
if ~isreal(values) || any(~isfinite(values)) || any(values < 0)
    error("radia:force:Geometry", "radius_m, axialLength_m, and angle_rad must be finite and nonnegative");
end
stress = radia.force.airGapShearStress( ...
    magneticFluxDensityRadial_T, magneticFluxDensityTangential_T, permeability_H_per_m);
torque_Nm = stress * radius_m^2 * angle_rad * axialLength_m;
end
