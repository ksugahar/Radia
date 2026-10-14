function torque_Nm = integrateTimeAverageMaxwellSurfaceTorque( ...
    magneticFluxDensityPhasor_T, outwardNormal, areaWeights_m2, ...
    samplePoints_m, pivot_m, permeability_H_per_m, amplitude)
%INTEGRATETIMEAVERAGEMAXWELLSURFACETORQUE Integrate phasor air-stress torque.

if nargin < 5
    pivot_m = [];
end
if nargin < 6
    permeability_H_per_m = [];
end
if nargin < 7
    amplitude = [];
end
[~, torque_Nm] = radia.force.integrateTimeAverageMaxwellSurfaceForceTorque( ...
    magneticFluxDensityPhasor_T, outwardNormal, areaWeights_m2, ...
    samplePoints_m, pivot_m, permeability_H_per_m, amplitude);
end
