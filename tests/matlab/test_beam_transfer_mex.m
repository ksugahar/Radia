function tests = test_beam_transfer_mex
%TEST_BEAM_TRANSFER_MEX Verify the native variational transfer-map ABI.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot,"matlab"));
testCase.TestData.SetupInfo = radia.setup(Force=true);
end

function testCommandCatalogAndNonlinearAttribution(testCase)
commands = string(radia.internal.callMex("api.commands"));
verifyTrue(testCase,ismember( ...
    "beam.transfer.propagate_variational",commands));

config = nonlinearConfig();
actual = runConfig(config);

verifyEqual(testCase,actual.schema, ...
    'radia.beam.variational-map.result.v1');
verifyEqual(testCase,actual.backend,'native-cpp-mex');
verifyEqual(testCase,string(actual.coordinate_order), ...
    ["x";"px_over_p0";"y";"py_over_p0";"sigma";"delta"]);
verifyEqual(testCase,actual.maximum_order,3);
verifySize(testCase,actual.R,[6,6]);
verifySize(testCase,actual.T,[6,6,6]);
verifySize(testCase,actual.U,[6,6,6,6]);
verifySize(testCase,actual.station_R,[6,6,4]);
verifySize(testCase,actual.region_T,[6,6,6,3]);
verifyEqual(testCase,actual.station_s_m,[0;0.3;0.7;0.9], ...
    "AbsTol",2e-15);
verifyEqual(testCase,string(actual.region_names), ...
    ["upstream_sextupole";"downstream_sextupole";"direct_octupole"]);
verifyEqual(testCase,actual.region_bounds_m, ...
    [0,0.3;0.3,0.7;0.7,0.9],"AbsTol",2e-15);

verifyEqual(testCase,actual.T(2,1,1),0.6,"AbsTol",2e-14);
verifyEqual(testCase,actual.T(3,1,2),0.6,"AbsTol",2e-14);
verifyEqual(testCase,actual.T(3,2,1),0.6,"AbsTol",2e-14);
verifyEqual(testCase,actual.U(3,1,1,1),1.08,"AbsTol",3e-14);
verifyEqual(testCase,actual.U(4,1,1,1),-0.14,"AbsTol",2e-14);
verifyEqual(testCase,actual.pair_regions,[1,2]);
verifyEqual(testCase,actual.pair_U_cascade(3,1,1,1,1),1.08, ...
    "AbsTol",3e-14);
verifyEqual(testCase,actual.region_U_direct(4,1,1,1,3),-0.14, ...
    "AbsTol",2e-14);

verifyLessThan(testCase,actual.diagnostics.R_composition_error,2e-14);
verifyLessThan(testCase,actual.diagnostics.T_reconstruction_error,3e-14);
verifyLessThan(testCase,actual.diagnostics.U_reconstruction_error,4e-14);
verifyLessThan(testCase,actual.diagnostics.T_symmetry_defect,1e-14);
verifyLessThan(testCase,actual.diagnostics.U_symmetry_defect,1e-14);
end

function testNormalQuadrupoleMatchesAnalyticMatrix(testCase)
strength = 1.7;
lengthM = 0.8;
rootStrength = sqrt(strength);
config.schema = 'radia.beam.variational-map.v1';
config.lengths_m = lengthM;
config.A_per_m = zeros(6,6);
config.A_per_m(1,2) = 1;
config.A_per_m(2,1) = -strength;
config.A_per_m(3,4) = 1;
config.A_per_m(4,3) = strength;
config.maximum_order = 1;
config.maximum_step_m = 0.001;

actual = radia.beam.propagateVariationalMap( ...
    config.lengths_m,config.A_per_m, ...
    MaximumOrder=config.maximum_order, ...
    MaximumStepM=config.maximum_step_m);
phase = rootStrength*lengthM;
expectedX = [cos(phase),sin(phase)/rootStrength; ...
    -rootStrength*sin(phase),cos(phase)];
expectedY = [cosh(phase),sinh(phase)/rootStrength; ...
    rootStrength*sinh(phase),cosh(phase)];
verifyEqual(testCase,actual.R(1:2,1:2),expectedX,"AbsTol",8e-13);
verifyEqual(testCase,actual.R(3:4,3:4),expectedY,"AbsTol",8e-13);
end

function testNonsymmetricJetAndInvalidShapeFailLoudly(testCase)
config.schema = 'radia.beam.variational-map.v1';
config.lengths_m = 1;
config.A_per_m = zeros(6,6);
config.F2_per_m = zeros(6,6,6);
config.F2_per_m(1,1,2) = 1;
config.maximum_order = 2;
verifyError(testCase,@() radia.internal.callMex( ...
    "beam.transfer.propagate_variational",config), ...
    "radia:mex:Exception");

config = rmfield(config,"F2_per_m");
config.A_per_m = zeros(5,6);
verifyError(testCase,@() radia.beam.propagateVariationalMap( ...
    config.lengths_m,config.A_per_m, ...
    MaximumOrder=config.maximum_order),"radia:beam:InvalidShape");
end

function testSubstepCascadeStaysInsideOneRegion(testCase)
f2 = zeros(6,6,6);
f2(2,1,1) = 2;
f2(3,1,2) = 1.5;
f2(3,2,1) = 1.5;
actual = radia.beam.propagateVariationalMap( ...
    1,zeros(6,6),F2PerM=f2,MaximumOrder=3,MaximumStepM=0.1);

verifyEqual(testCase,actual.U(3,1,1,1),4.5,"AbsTol",2e-13);
verifyEqual(testCase,actual.region_U_local_cascade(3,1,1,1,1), ...
    4.5,"AbsTol",2e-13);
verifyTrue(testCase,isempty(actual.pair_regions));
verifyTrue(testCase,isempty(actual.pair_U_cascade));
verifyLessThan(testCase,actual.diagnostics.U_reconstruction_error,3e-13);
end

function testUnicodeRegionNameRoundTrips(testCase)
regionName = string(char([20837,21475,22235,26997]));
actual = radia.beam.propagateVariationalMap( ...
    0.1,zeros(6),Names=regionName,MaximumOrder=1,MaximumStepM=0.1);
verifyEqual(testCase,string(actual.region_names{1}),regionName);
end

function config = nonlinearConfig()
config.schema = 'radia.beam.variational-map.v1';
config.lengths_m = [0.3,0.4,0.2];
config.A_per_m = zeros(6,6,3);
config.F2_per_m = zeros(6,6,6,3);
config.F3_per_m = zeros(6,6,6,6,3);
config.F2_per_m(2,1,1,1) = 2;
config.F2_per_m(3,1,2,2) = 1.5;
config.F2_per_m(3,2,1,2) = 1.5;
config.F3_per_m(4,1,1,1,3) = -0.7;
config.names = { ...
    'upstream_sextupole', ...
    'downstream_sextupole', ...
    'direct_octupole'};
config.maximum_order = 3;
config.maximum_step_m = 1;
end

function result = runConfig(config)
result = radia.beam.propagateVariationalMap( ...
    config.lengths_m,config.A_per_m, ...
    F2PerM=config.F2_per_m,F3PerM=config.F3_per_m, ...
    Names=config.names,MaximumOrder=config.maximum_order, ...
    MaximumStepM=config.maximum_step_m);
end
