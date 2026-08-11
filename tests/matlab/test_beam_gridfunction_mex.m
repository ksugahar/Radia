function tests = test_beam_gridfunction_mex
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabRoot = fullfile(repoRoot,"matlab");
if ~contains(path,matlabRoot)
    addpath(matlabRoot);
end
testCase.TestData.VolPath = fullfile( ...
    repoRoot,"tests","fixtures","beam","affine_field_tetra.vol");
end

function testCommandIsPublished(testCase)
commands = string(radia.internal.callMex("api.commands"));
verifyTrue(testCase,any(commands == "beam.transfer.from_grid_function"));
end

function testCoordinateGridFunctionFeedsNativeLinearMap(testCase)
field = radia.ngsolve.GridFunction.create( ...
    string(testCase.TestData.VolPath),"vectorh1",1);
fieldCleanup = onCleanup(@() delete(field));
coordinates = radia.ngsolve.CoefficientFunction.coordinates(3);
coordinateCleanup = onCleanup(@() delete(coordinates));
field.interpolate(coordinates);

lengths = [0.2;0.3];
positions = [0 0.2 0.25;0 0.2 0.75];
tangents = repmat([0 0 1],2,1);
result = radia.beam.propagateGridFunctionLinearMap( ...
    field,lengths,positions,tangents,2, ...
    SampleRadiusM=0.02,Names=["entrance";"出口四極"]);

verifyEqual(testCase,string(result.schema), ...
    "radia.beam.grid-function-linear-map.result.v1");
verifyEqual(testCase,string(result.backend), ...
    "native-cpp-ngsolve-gridfunction-mex");
verifyEqual(testCase,result.maximum_order,1);
verifyEqual(testCase,result.curvature_per_m,0.1*ones(2,1), ...
    "AbsTol",1e-12);
verifyEqual(testCase,result.normal_gradient_per_m2,zeros(2,1), ...
    "AbsTol",1e-12);
verifyEqual(testCase,result.skew_gradient_per_m2,zeros(2,1), ...
    "AbsTol",1e-12);
verifyEqual(testCase,result.transverse_divergence_t_per_m,2*ones(2,1), ...
    "AbsTol",1e-11);
verifyEqual(testCase,result.fit_rank,3*ones(2,1));
verifyEqual(testCase,result.scaled_design_condition,1.5*ones(2,1));
verifyEqual(testCase,result.center_fit_bias_t,zeros(2,1),"AbsTol",1e-12);
verifyLessThan(testCase,max(result.maximum_fit_residual_t),1e-12);
verifyEqual(testCase,string(result.region_names(:)), ...
    ["entrance";"出口四極"]);

generator = zeros(6);
generator(1,2) = 1;
generator(2,1) = -0.1^2;
generator(2,6) = 0.1;
generator(3,4) = 1;
generator(5,1) = 0.1;
verifyEqual(testCase,result.local_A_per_m, ...
    repmat(generator,1,1,2),"AbsTol",1e-12);
verifyEqual(testCase,result.R,expm(generator*sum(lengths)), ...
    "AbsTol",2e-11);

clear fieldCleanup coordinateCleanup
end

function testOutsideMeshFailsLoudly(testCase)
field = radia.ngsolve.GridFunction.create( ...
    string(testCase.TestData.VolPath),"vectorh1",1);
fieldCleanup = onCleanup(@() delete(field));
coordinates = radia.ngsolve.CoefficientFunction.coordinates(3);
coordinateCleanup = onCleanup(@() delete(coordinates));
field.interpolate(coordinates);

verifyError(testCase,@() radia.beam.propagateGridFunctionLinearMap( ...
    field,0.1,[0 0 4],[0 0 1],2),"radia:mex:Exception");

clear fieldCleanup coordinateCleanup
end
