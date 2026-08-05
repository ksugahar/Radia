function tests = test_maglev_model
%TEST_MAGLEV_MODEL Verify the standalone moving HCurl/CLN MagLev model.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.FileGenConfig = Simulink.fileGenControl("getConfig");
testCase.TestData.FileGenRoot = string(tempname("C:\temp"));
Simulink.fileGenControl("set", ...
    CacheFolder=fullfile(testCase.TestData.FileGenRoot, "cache"), ...
    CodeGenFolder=fullfile(testCase.TestData.FileGenRoot, "codegen"), ...
    createDir=true);
end

function teardownOnce(testCase)
Simulink.fileGenControl("setConfig", ...
    config=testCase.TestData.FileGenConfig);
if isfolder(testCase.TestData.FileGenRoot)
    rmdir(testCase.TestData.FileGenRoot, "s");
end
end

function testSmokeFamilyIsPassiveAndPositionDependent(testCase)
family = radia.simulink.makeMagLevSmokeFamily(SampleTime_s=1.0e-3);
verifyEqual(testCase, family.schema, "radia.hcurl.eddy_cln.family.v1");
verifyTrue(testCase, family.shared_state_basis);
verifyEqual(testCase, family.snapshot_count, 2);
verifyEqual(testCase, family.sample_time_s, 1.0e-3, "AbsTol", 0);
verifyNotEqual(testCase, family.models{1}.Ad, family.models{2}.Ad);
verifyTrue(testCase, all(cellfun(@(model) model.passive, family.models)));
end

function testBuilderCreatesRunnableModel(testCase)
outputDirectory = string(tempname("C:\temp"));
mkdir(outputDirectory);
modelName = "radia_maglev_model_" + ...
    erase(string(java.util.UUID.randomUUID), "-");
cleanup = onCleanup(@() closeIfLoaded(modelName));
modelPath = radia.simulink.buildMagLevModel( ...
    ModelName=modelName, OutputDirectory=outputDirectory, Open=false);
verifyTrue(testCase, isfile(modelPath));
load_system(modelPath);

plant = modelName + "/MagLev Plant";
verifyEqual(testCase, string(get_param(plant, "Mask")), "on");
plantMask = Simulink.Mask.get(plant);
familyParameter = plantMask.getParameter("family");
verifyEqual(testCase, string(familyParameter.Evaluate), "off");
verifyTrue(testCase, contains(string(familyParameter.Callback), ...
    "onMagLevBlockFamilyChanged"));
familyExpression = ...
    "radia.simulink.makeMagLevSmokeFamily(SampleTime_s=1e-4)";
set_param(plant, "family", familyExpression);
radia.simulink.onMagLevBlockFamilyChanged(plant);
verifyEqual(testCase, string(get_param( ...
    plant + "/Moving HCurl CLN", "Parameters")), familyExpression);
verifyEqual(testCase, string(get_param( ...
    plant + "/Moving HCurl CLN", "FunctionName")), ...
    "radia_hcurl_eddy_cln_family_sfunction");
ports = get_param(plant, "PortHandles");
verifyEqual(testCase, numel(ports.Inport), 3);
verifyEqual(testCase, numel(ports.Outport), 2);
contract = get_param(plant, "UserData");
verifyFalse(testCase, contract.python_per_step);
verifyFalse(testCase, contract.surrogate);
verifyEqual(testCase, string(get_param(modelName, "Solver")), ...
    "FixedStepDiscrete");

set_param(modelName, "SimulationCommand", "update");
simulation = sim(modelName, "StopTime", "0.01", ...
    "ReturnWorkspaceOutputs", "on");
outputs = simulation.get("yout");
verifyEqual(testCase, outputs.numElements, 2);
force = outputs.getElement(2).Values.Data;
verifySize(testCase, force, [101, 3]);
verifyTrue(testCase, all(isfinite(force), "all"));
verifyGreaterThan(testCase, max(abs(force(:,3))), 0);
clear cleanup
closeIfLoaded(modelName);
end

function testParameterMaskResamplesEmbeddedFamily(testCase)
modelName = "radia_maglev_callback_" + ...
    erase(string(java.util.UUID.randomUUID), "-");
outputDirectory = string(tempname("C:\temp"));
mkdir(outputDirectory);
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildMagLevModel( ...
    ModelName=modelName, OutputDirectory=outputDirectory, Open=false);
load_system(fullfile(outputDirectory, modelName + ".slx"));
parameterPath = modelName + "/MagLev Parameters";
set_param(parameterPath, "sample_time_s", "0.002", ...
    "interpolation", "nearest", "extrapolation", "clamp");
radia.simulink.onMagLevFamilyChanged(parameterPath);
workspace = get_param(modelName, "ModelWorkspace");
family = workspace.getVariable("radia_maglev_family");
verifyEqual(testCase, family.sample_time_s, 0.002, "AbsTol", 0);
verifyEqual(testCase, string(family.interpolation), "nearest");
verifyEqual(testCase, string(family.extrapolation), "clamp");
verifyEqual(testCase, string(get_param(modelName, "FixedStep")), "0.002");
clear cleanup
closeIfLoaded(modelName);
end

function testTrackedModelLoadsAndUpdates(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
modelPath = fullfile(root, "matlab", "radia_maglev.slx");
verifyTrue(testCase, isfile(modelPath));
load_system(modelPath);
cleanup = onCleanup(@() closeIfLoaded("radia_maglev"));
set_param("radia_maglev", "SimulationCommand", "update");
verifyEqual(testCase, string(get_param( ...
    "radia_maglev/MagLev Plant/Moving HCurl CLN", "FunctionName")), ...
    "radia_hcurl_eddy_cln_family_sfunction");
verifyEqual(testCase, string(get_param( ...
    "radia_maglev/MagLev Parameters", "Mask")), "on");
clear cleanup
closeIfLoaded("radia_maglev");
end

function closeIfLoaded(name)
if bdIsLoaded(name)
    close_system(name, 0);
end
end
