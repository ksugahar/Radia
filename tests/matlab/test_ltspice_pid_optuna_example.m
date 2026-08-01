function tests = test_ltspice_pid_optuna_example
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
testCase.TestData.Root = root;
end

function teardownOnce(testCase)
rmpath(fullfile(testCase.TestData.Root,"matlab"));
end

function testGeneratedModelsUseRadiaBlocks(testCase)
root = testCase.TestData.Root;
files = struct("plant",fullfile(root,"matlab","radia_ltspice_pid_plant.slx"), ...
    "harness",fullfile(root,"matlab","radia_ltspice_pid_optuna.slx"));
cleanup = onCleanup(@() closeModels());
load_system(files.plant); load_system(files.harness);
verifyEqual(testCase,string(get_param( ...
    "radia_ltspice_pid_plant/LTspice Circuit","FunctionName")), ...
    "radia_ltspice_sfun");
verifyEqual(testCase,string(get_param( ...
    "radia_ltspice_pid_optuna/Optuna Optimization","FunctionName")), ...
    "radia_optuna_sfun");
verifyEqual(testCase,string(get_param( ...
    "radia_ltspice_pid_optuna/Optuna Optimization","objective_fcn")), ...
    "'radia.simulink.ltspicePIDObjective'");
clear cleanup; closeModels();
end

function closeModels
models = ["radia_ltspice_pid_plant","radia_ltspice_pid_optuna"];
for model = models
    if bdIsLoaded(model), close_system(model,0); end
end
end
