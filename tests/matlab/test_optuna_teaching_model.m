function tests = test_optuna_teaching_model
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root, "matlab"));
testCase.TestData.Root = root;
end

function testAllExercisesBuildAndUpdateCleanly(testCase)
exercises = ["quadratic","pareto","reliability"];
for exercise = exercises
    folder = string(tempname("C:\temp"));
    mkdir(folder);
    cleanup = onCleanup(@()rmdir(folder, "s"));
    path = folder + "/teaching_" + exercise + ".slx";
    actual = radia.simulink.buildOptunaTeachingModel( ...
        OutputPath=path, Exercise=exercise, ...
        StoragePath=folder + "/study.mat");
    verifyEqual(testCase, replace(actual, "/", "\"), ...
        replace(path, "/", "\"));
    verifyTrue(testCase, isfile(path));
    [~, name] = fileparts(path);
    load_system(path);
    modelCleanup = onCleanup(@()close_system(name, 0));
    verifyEqual(testCase, replace(string(get_param(name, "FileName")), ...
        "/", "\"), replace(path, "/", "\"));
    set_param(name, "SimulationCommand", "update");
    ports = get_param(name + "/Optuna Study", "PortHandles");
    verifyEqual(testCase, numel(ports.Inport), 2);
    verifyEqual(testCase, numel(ports.Outport), 4);
    runtimePorts = get_param( ...
        name + "/Optuna Study/Advanced Runtime", "PortHandles");
    verifyEqual(testCase, numel(runtimePorts.Inport), 6);
    verifyEqual(testCase, numel(runtimePorts.Outport), 18);
    contract = get_param(name + "/Optuna Study", "UserData");
    verifyEqual(testCase, string(contract.interface), "compact-student");
    verifyEqual(testCase, string(contract.iteration), ...
        "mask-edit;resume-study;review-table;apply-trial");
    workspace = get_param(name, "ModelWorkspace");
    verifyTrue(testCase, workspace.hasVariable("x"));
    verifyTrue(testCase, contains(string(get_param( ...
        name + "/Optuna Study", "parameter_spec")), ...
        "radia.optuna.OptimizationParameter"));
    clear modelCleanup cleanup
end
end

function testQuadraticExerciseRunsToCompletion(testCase)
[out, review, cleanup] = runExercise("quadratic"); %#ok<ASGLU>
status = out.get("teaching_status");
attempted = out.get("teaching_attempted");
best = out.get("teaching_best");
verifyEqual(testCase, status(end), 1);
verifyEqual(testCase, attempted(end), 12);
verifyTrue(testCase, isfinite(best(end)));
verifyEqual(testCase, height(review.trials), 12);
clear cleanup
end

function testParetoExercisePublishesFront(testCase)
[out, review, cleanup] = runExercise("pareto"); %#ok<ASGLU>
status = out.get("teaching_status");
verifyEqual(testCase, status(end), 1);
verifyGreaterThan(testCase, height(review.pareto), 1);
clear cleanup
end

function testReliabilityExerciseRecordsPrunedAndFailed(testCase)
[out, review, cleanup] = runExercise("reliability"); %#ok<ASGLU>
status = out.get("teaching_status");
attempted = out.get("teaching_attempted");
verifyEqual(testCase, status(end), 1);
verifyEqual(testCase, attempted(end), 5);
verifyGreaterThan(testCase, sum(review.trials.State == "PRUNED"), 0);
verifyGreaterThan(testCase, sum(review.trials.State == "FAIL"), 0);
clear cleanup
end

function testMaskIterationNeedsNoRewiring(testCase)
folder = string(tempname("C:\temp"));
mkdir(folder);
modelPath = folder + "/iteration.slx";
studyA = folder + "/study-a.mat";
studyB = folder + "/study-b.mat";
radia.simulink.buildOptunaTeachingModel( ...
    OutputPath=modelPath, Exercise="quadratic", StoragePath=studyA);
load_system(modelPath);
cleanup = onCleanup(@()cleanExercise("iteration", folder, ""));
block = "iteration/Optuna Study";
portsBefore = get_param(block, "PortConnectivity");

sim("iteration", ReturnWorkspaceOutputs="on");
first = radia.simulink.reviewOptunaStudy(block);
verifyEqual(testCase, height(first.trials), 12);

comparisonSpace = ...
    "radia.optuna.OptimizationParameter('x',Value=0," + ...
    "Minimum=-0.5,Maximum=0.5)";
set_param(block, "parameter_spec", comparisonSpace);
verifyError(testCase, ...
    @()sim("iteration", ReturnWorkspaceOutputs="on"), ...
    "Simulink:blocks:MSFB_BlockMethodFailed");
unchanged = radia.simulink.reviewOptunaStudy(block);
verifyEqual(testCase, height(unchanged.trials), 12);

set_param(block, "num_trials", "8", "sampler_name", "qmc", ...
    "seed", "17", "pruner_name", "none", ...
    "storage_path", "'" + studyB + "'", ...
    "parameter_spec", comparisonSpace);
sim("iteration", ReturnWorkspaceOutputs="on");
comparison = radia.simulink.reviewOptunaStudy(block);
verifyEqual(testCase, height(comparison.trials), 8);

set_param(block, "num_trials", "12");
out = sim("iteration", ReturnWorkspaceOutputs="on");
resumed = radia.simulink.reviewOptunaStudy(block);
verifyEqual(testCase, height(resumed.trials), 12);
verifyEqual(testCase, out.teaching_status(end), 1);

values = radia.simulink.applyOptunaTrial(block, -1);
workspace = get_param("iteration", "ModelWorkspace");
verifyEqual(testCase, getVariable(workspace, "x"), values.x, ...
    "AbsTol", 0);
portsAfter = get_param(block, "PortConnectivity");
verifyEqual(testCase, portsAfter, portsBefore);
clear cleanup
end

function [out, review, cleanup] = runExercise(exercise)
folder = string(tempname("C:\temp"));
mkdir(folder);
path = folder + "/run_" + exercise + ".slx";
radia.simulink.buildOptunaTeachingModel( ...
    OutputPath=path, Exercise=exercise, ...
    StoragePath=folder + "/study.mat");
[~, name] = fileparts(path);
cleanup = onCleanup(@()cleanExercise(name, folder, exercise));
load_system(path);
out = sim(name, ReturnWorkspaceOutputs="on");
review = radia.simulink.reviewOptunaStudy(name + "/Optuna Study");
end

function cleanExercise(name, folder, ~)
if bdIsLoaded(name)
    close_system(name, 0);
end
if isfolder(folder)
    rmdir(folder, "s");
end
end
