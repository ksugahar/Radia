function torque_Nm = integrateMaxwellSurfaceTorque( ...
    magneticFluxDensity_T, outwardNormal, areaWeights_m2, samplePoints_m, ...
    pivot_m, permeability_H_per_m)
%INTEGRATEMAXWELLSURFACETORQUE Integrate static air-stress torque.

if nargin < 5
    pivot_m = [];
end
if nargin < 6
    permeability_H_per_m = [];
end
[~, torque_Nm] = radia.force.integrateMaxwellSurfaceForceTorque( ...
    magneticFluxDensity_T, outwardNormal, areaWeights_m2, samplePoints_m, ...
    pivot_m, permeability_H_per_m);
end
