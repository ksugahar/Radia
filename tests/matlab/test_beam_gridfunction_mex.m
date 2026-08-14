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
expected = ["beam.transfer.from_grid_function"; ...
    "beam.transfer.multipole_from_grid_function"; ...
    "beam.track.grid_function"];
verifyTrue(testCase,all(ismember(expected,commands)));
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

function testMultipoleWrapperPublishesRequestedFitAndMapOrder(testCase)
field = radia.ngsolve.GridFunction.create( ...
    string(testCase.TestData.VolPath),"vectorh1",1);
fieldCleanup = onCleanup(@() delete(field));
constantField = radia.ngsolve.CoefficientFunction.constant([0;0.2;0]);
constantCleanup = onCleanup(@() delete(constantField));
field.interpolate(constantField);

result = radia.beam.propagateGridFunctionMultipoleMap( ...
    field,0.02,[0 0.2 0.25],[0 0 1],2, ...
    SampleRadiusM=0.02,MultipoleOrder=2,MaximumMapOrder=1);
verifyEqual(testCase,string(result.schema), ...
    "radia.beam.grid-function-multipole-map.result.v1");
verifyEqual(testCase,result.maximum_order,1);
verifyEqual(testCase,result.linearization_order,2);
verifyEqual(testCase,string(result.fit_model), ...
    "nine-point transverse harmonic multipole expansion through order 2");
verifyEqual(testCase,result.multipole_normal_t_per_m_power(1,1),0.2, ...
    "AbsTol",1e-12);

nonlinear = radia.beam.propagateGridFunctionMultipoleMap( ...
    field,0.02,[0 0.2 0.25],[0 0 1],2, ...
    SampleRadiusM=0.02,MultipoleOrder=3,MaximumMapOrder=3, ...
    MaximumStepM=1e-3);
verifyEqual(testCase,nonlinear.maximum_order,3);
verifyGreaterThan(testCase,max(abs(nonlinear.T(:))),0);
verifyGreaterThan(testCase,max(abs(nonlinear.U(:))),0);

clear fieldCleanup constantCleanup
end

function testGridFunctionTrackingUsesNativePointEvaluation(testCase)
field = radia.ngsolve.GridFunction.create( ...
    string(testCase.TestData.VolPath),"vectorh1",1);
fieldCleanup = onCleanup(@() delete(field));
constantField = radia.ngsolve.CoefficientFunction.constant([0;0.2;0]);
constantCleanup = onCleanup(@() delete(constantField));
field.interpolate(constantField);

species = radia.beam.proton();
reference = radia.beam.referenceParticle(species,220e6);
state = radia.beam.cartesianState( ...
    [0 0.2 0.25],[0 0 reference.momentum_kg_m_s]);
trajectory = radia.beam.trackGridFunction( ...
    species,state,field,0,0.005,5e-4);

verifyEqual(testCase,string(trajectory.backend), ...
    "native-cpp-ngsolve-gridfunction-mex");
verifyGreaterThan(testCase,trajectory.summary.accepted_steps,0);
finalState = trajectory.samples{end};
verifyEqual(testCase,finalState.path_length_m,0.005,"AbsTol",1e-14);
verifyLessThan(testCase,finalState.position_m(1),0);

clear fieldCleanup constantCleanup
end
