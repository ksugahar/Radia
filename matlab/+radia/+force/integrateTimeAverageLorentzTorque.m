function torque_Nm = integrateTimeAverageLorentzTorque( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, ...
    volumeWeights_m3, samplePoints_m, pivot_m, amplitude)
%INTEGRATETIMEAVERAGELORENTZTORQUE Integrate phasor Lorentz torque.

if nargin < 5
    pivot_m = [];
end
if nargin < 6
    amplitude = [];
end
[~, torque_Nm] = radia.force.integrateTimeAverageLorentzForceTorque( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, ...
    volumeWeights_m3, samplePoints_m, pivot_m, amplitude);
end
