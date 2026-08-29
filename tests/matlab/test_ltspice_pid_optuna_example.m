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
    "radia_ltspice_pid_optuna/Optuna Study/Advanced Runtime", ...
    "FunctionName")), ...
    "radia_optuna_sfun");
verifyEqual(testCase,string(get_param( ...
    "radia_ltspice_pid_optuna/Optuna Study","objective_fcn")), ...
    "'radia.simulink.ltspicePIDObjective'");
clear cleanup; closeModels();
end

function testGeneratorOwnsPlantAndUsesBundledNetlist(testCase)
outputDirectory=fullfile("C:\temp","radia_ltspice_pid_generator_test");
files=radia.simulink.buildLTspicePIDOptunaExample( ...
    OutputDirectory=outputDirectory,NumTrials=2,OpenModel=false);
cleanup=onCleanup(@() closeModels());
verifyTrue(testCase,isfile(files.plant));
verifyTrue(testCase,isfile(files.harness));
load_system(files.plant);load_system(files.harness);
netlist=fullfile(testCase.TestData.Root,"matlab","samples", ...
    "ltspice_pid_rc_plant.cir");
verifyTrue(testCase,isfile(netlist));
verifySubstring(testCase,string(get_param( ...
    "radia_ltspice_pid_plant/LTspice Circuit","Parameters")),netlist);
objective=string(get_param( ...
    "radia_ltspice_pid_optuna/Optuna Study","objective_fcn"));
verifyEqual(testCase,objective,"'radia.simulink.ltspicePIDObjective'");
verifyFalse(testCase,contains(objective,string(outputDirectory)));
clear cleanup;closeModels();
end

function closeModels
models = ["radia_ltspice_pid_plant","radia_ltspice_pid_optuna"];
for model = models
    if bdIsLoaded(model), close_system(model,0); end
end
end
