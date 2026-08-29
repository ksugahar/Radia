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
    ports = get_param(name + "/Optuna Optimization", "PortHandles");
    verifyEqual(testCase, numel(ports.Inport), 6);
    verifyEqual(testCase, numel(ports.Outport), 18);
    workspace = get_param(name, "ModelWorkspace");
    parameters = getVariable(workspace, ...
        "radia_optuna_teaching_parameters");
    verifyClass(testCase, parameters, ...
        "radia.optuna.OptimizationParameter");
    clear modelCleanup cleanup
end
end

function testQuadraticExerciseRunsToCompletion(testCase)
[out, cleanup] = runExercise("quadratic"); %#ok<ASGLU>
status = out.get("teaching_status");
attempted = out.get("teaching_attempted");
checkpoint = out.get("teaching_checkpoint");
best = out.get("teaching_best");
verifyEqual(testCase, status(end), 1);
verifyEqual(testCase, attempted(end), 12);
verifyGreaterThan(testCase, checkpoint(end), 0);
verifyTrue(testCase, isfinite(best(end)));
clear cleanup
end

function testParetoExercisePublishesFront(testCase)
[out, cleanup] = runExercise("pareto"); %#ok<ASGLU>
status = out.get("teaching_status");
paretoCount = out.get("teaching_pareto_count");
verifyEqual(testCase, status(end), 1);
verifyGreaterThan(testCase, paretoCount(end), 1);
clear cleanup
end

function testReliabilityExerciseRecordsPrunedAndFailed(testCase)
[out, cleanup] = runExercise("reliability"); %#ok<ASGLU>
status = out.get("teaching_status");
attempted = out.get("teaching_attempted");
pruned = out.get("teaching_pruned");
failed = out.get("teaching_failed");
verifyEqual(testCase, status(end), 1);
verifyEqual(testCase, attempted(end), 5);
verifyGreaterThan(testCase, pruned(end), 0);
verifyGreaterThan(testCase, failed(end), 0);
clear cleanup
end

function [out, cleanup] = runExercise(exercise)
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
end

function cleanExercise(name, folder, ~)
if bdIsLoaded(name)
    close_system(name, 0);
end
if isfolder(folder)
    rmdir(folder, "s");
end
end
