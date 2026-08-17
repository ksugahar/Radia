function tests = test_electromagnetic_force
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(root, "matlab");
entries = string(strsplit(path, pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries, string(matlabDirectory)));
testCase.TestData.Path = matlabDirectory;
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.Path);
end
end

function testLorentzDensityAndIntegral(testCase)
density = radia.force.lorentzDensity( ...
    [0, 0, 2e6], [0, 0.3, 0]);
verifyEqual(testCase, density, [-6e5, 0, 0], "AbsTol", 1e-12);
force = radia.force.integrateLorentz( ...
    [0, 0, 2e6], [0, 0.3, 0], 2.5e-6);
verifyEqual(testCase, force, [-1.5, 0, 0], "AbsTol", 1e-14);
end

function testNormalFieldMaxwellPressure(testCase)
mu0 = 4*pi*1e-7;
pressure = 1.2^2/(2*mu0);
stress = radia.force.maxwellStressAir([0, 0, 1.2]);
traction = radia.force.maxwellTractionAir( ...
    [0, 0, 1.2], [0, 0, 4]);
verifyEqual(testCase, stress, diag([-pressure, -pressure, pressure]), ...
    "RelTol", 2e-15);
verifyEqual(testCase, traction, [0, 0, pressure], "RelTol", 2e-15);
end

function testTwoWireLorentzAndMaxwellRoutes(testCase)
mu0 = 4*pi*1e-7;
sourceCurrent = 80;
targetCurrent = 35;
separation = 0.06;
targetRadius = 0.008;
expected = -mu0*sourceCurrent*targetCurrent/(2*pi*separation);
targetArea = pi*targetRadius^2;
lorentz = radia.force.integrateLorentz( ...
    [0, 0, targetCurrent/targetArea], ...
    [0, mu0*sourceCurrent/(2*pi*separation), 0], targetArea);

sampleCount = 4096;
angle = 2*pi*((0:sampleCount-1).' + 0.5)/sampleCount;
normals = [cos(angle), sin(angle), zeros(sampleCount, 1)];
targetCenter = [separation, 0, 0];
points = targetCenter + targetRadius*normals;
field = localLineCurrentField(points, [0, 0, 0], sourceCurrent) ...
    + localLineCurrentField(points, targetCenter, targetCurrent);
weights = repmat(targetRadius*2*pi/sampleCount, sampleCount, 1);
maxwell = radia.force.integrateMaxwellSurface(field, normals, weights);

verifyEqual(testCase, lorentz, [expected, 0, 0], ...
    "RelTol", 1e-14, "AbsTol", 1e-18);
verifyEqual(testCase, maxwell, [expected, 0, 0], ...
    "RelTol", 2e-13, "AbsTol", 1e-18);
end

function testStaticForceTorqueResultants(testCase)
[force, torque] = radia.force.integrateLorentzForceTorque( ...
    [0, 0, 1], [0, 1, 0], 1, [0, 2, 0]);
verifyEqual(testCase, force, [-1, 0, 0], "AbsTol", 1e-14);
verifyEqual(testCase, torque, [0, 0, 2], "AbsTol", 1e-14);
verifyEqual(testCase, radia.force.integrateLorentzTorque( ...
    [0, 0, 1], [0, 1, 0], 1, [0, 2, 0]), torque, "AbsTol", 1e-14);

field = [1, 0, 0];
normal = [1, 0, 0];
[surfaceForce, surfaceTorque] = radia.force.integrateMaxwellSurfaceForceTorque( ...
    field, normal, 2, [0, 1, 0]);
verifyEqual(testCase, surfaceTorque, cross([0, 1, 0], surfaceForce), ...
    "RelTol", 2e-15);
end

function testPeakAndRmsPhasorConventions(testCase)
staticDensity = radia.force.lorentzDensity([0, 0, 4], [0, 3, 0]);
peakDensity = radia.force.timeAverageLorentzDensity( ...
    [0, 0, 4], [0, 3, 0], "peak");
rmsDensity = radia.force.timeAverageLorentzDensity( ...
    [0, 0, 4], [0, 3, 0], "rms");
verifyEqual(testCase, peakDensity, 0.5*staticDensity, "AbsTol", 1e-14);
verifyEqual(testCase, rmsDensity, staticDensity, "AbsTol", 1e-14);
verifyEqual(testCase, radia.force.timeAverageLorentzDensity( ...
    [0, 0, 4i], [0, 3, 0], "peak"), [0, 0, 0], "AbsTol", 1e-14);

staticStress = radia.force.maxwellStressAir([0, 0, 1.2]);
peakStress = radia.force.timeAverageMaxwellStressAir( ...
    [0, 0, 1.2], [], "peak");
verifyEqual(testCase, peakStress, 0.5*staticStress, "RelTol", 2e-15);

[force, torque] = radia.force.integrateTimeAverageLorentzForceTorque( ...
    [0, 0, 4], [0, 3, 0], 1, [0, 2, 0], [], "peak");
verifyEqual(testCase, force, 0.5*staticDensity, "AbsTol", 1e-14);
verifyEqual(testCase, torque, cross([0, 2, 0], force), "AbsTol", 1e-14);
end

function testVirtualWorkCoenergyAndAirGapTorque(testCase)
positions = linspace(-0.002, 0.002, 5).';
forceExpected = 7.5;
verifyEqual(testCase, radia.force.virtualWorkForce( ...
    positions, 0.25 + forceExpected*positions), ...
    repmat(forceExpected, 5, 1), "AbsTol", 1e-12);

angles = linspace(0, 0.4, 5).';
torqueExpected = -2.25;
verifyEqual(testCase, radia.force.coenergyTorque( ...
    angles, 1 + torqueExpected*angles), ...
    repmat(torqueExpected, 5, 1), "AbsTol", 1e-12);
verifyEqual(testCase, radia.force.airGapShearTorque( ...
    0.8, 0.1, 0.05, 0.1), 100, "RelTol", 2e-15);

angles = (0:3).' * (pi/2);
static = radia.force.airGapShearTorqueFromAngleSamples( ...
    angles, repmat(0.8, 4, 1), repmat(0.1, 4, 1), 0.05, 0.1);
phasor = radia.force.timeAverageAirGapShearTorqueFromAngleSamples( ...
    angles, complex(repmat(0.8, 4, 1)), complex(repmat(0.1, 4, 1)), ...
    0.05, 0.1, true, 2*pi, [], "peak");
verifyEqual(testCase, static.integrated_angle_rad, 2*pi, "RelTol", 2e-15);
verifyEqual(testCase, static.torque_Nm, 100, "RelTol", 2e-15);
verifyEqual(testCase, phasor.torque_Nm, 50, "RelTol", 2e-15);
end

function testForceTorqueResultContract(testCase)
result = radia.force.forceTorqueResult( ...
    [1, 2, 3], [0.1, 0.2, 0.3], ...
    Method="time_average_lorentz", Frame="rotor", ...
    Pivot_m=[0, 0, 0.1], FieldConvention="time_average_phasor", ...
    Amplitude="peak");
verifyEqual(testCase, result.schema, "radia.force-result/v1");
verifyEqual(testCase, result.frame, "rotor");
verifyEqual(testCase, result.force_N, [1, 2, 3]);
verifyEqual(testCase, result.torque_Nm, [0.1, 0.2, 0.3]);
verifyEqual(testCase, result.phasor_amplitude, "peak");
end

function testInvalidContractsFailLoudly(testCase)
verifyError(testCase, @()radia.force.lorentzDensity( ...
    [1, 0], [0, 0, 1]), "radia:force:VectorShape");
verifyError(testCase, @()radia.force.maxwellTractionAir( ...
    [0, 0, 1], [0, 0, 0]), "radia:force:Normal");
verifyError(testCase, @()radia.force.integrateLorentz( ...
    [0, 0, 1; 0, 0, 1], [0, 1, 0], [1; -1]), ...
    "radia:force:Weights");
end

function field = localLineCurrentField(points, center, current)
relative = points - center;
radiusSquared = relative(:, 1).^2 + relative(:, 2).^2;
scale = (4*pi*1e-7)*current./(2*pi*radiusSquared);
field = [-scale.*relative(:, 2), scale.*relative(:, 1), ...
    zeros(size(points, 1), 1)];
end
