function tests = test_maglev_model
%TEST_MAGLEV_MODEL Verify the standalone moving Foster HCurl MagLev model.
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

function testSmokeFamilySharesModesAndMovesPorts(testCase)
family = radia.simulink.makeMagLevSmokeFamily(SampleTime_s=1.0e-3);
verifyEqual(testCase, family.schema, "radia.hcurl.eddy_foster.family.v1");
verifyTrue(testCase, family.shared_modes);
verifyEqual(testCase, family.snapshot_count, 2);
verifyEqual(testCase, family.sample_time_s, 1.0e-3, "AbsTol", 0);
verifyGreaterThan(testCase, min(family.decay_rates), 0);
verifyEqual(testCase, diag(family.Ad), exp(-family.decay_rates * 1.0e-3), ...
    "RelTol", 1e-15);
verifyNotEqual(testCase, family.modal_port_rhs(:,:,1), family.modal_port_rhs(:,:,2));
middle = radia.simulink.interpolateHCurlEddyFosterFamily(family, 0.005);
verifyEqual(testCase, middle.modal_port_rhs, ...
    0.5 * (family.modal_port_rhs(:,:,1) + family.modal_port_rhs(:,:,2)), "AbsTol", 1e-15);
verifyEqual(testCase, middle.Ad, family.Ad, "AbsTol", 0);
verifyError(testCase, @() radia.simulink.interpolateHCurlEddyFosterFamily( ...
    family, 0.02), "radia:simulink:HCurlFosterExtrapolation");
end

function testTimeDomainForceAveragesToHarmonicForce(testCase)
% Instantaneous block force K*z(t)*i(t) must average to the phasor force.
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
familyFile = fullfile(root, "validation_test", "maglev", ...
    "team28_coilbuilder_hcurl_eddy_foster_family.json");
sampleTime = 1.0e-5;
family = radia.simulink.loadHCurlEddyFosterFamily(familyFile, SampleTime_s=sampleTime);
model = radia.simulink.interpolateHCurlEddyFosterFamily(family, 0.0);
omega = 2 * pi * 50;
amplitude = 20;
steps = round(0.1 / sampleTime);
z = zeros(family.state_order, 1);
lift = zeros(steps, 1);
for k = 1:steps
    t = (k - 1) * sampleTime;
    force = radia.internal.hcurlEddyFosterInstantaneousForce( ...
        model.force_operator, z, amplitude * sin(omega * t));
    lift(k) = force(3);
    z = model.Ad * z + model.Bd * (-amplitude * omega * cos(omega * t));
end
period = round(0.02 / sampleTime);
harmonic = radia.simulink.evaluateHCurlEddyFosterForce(model, ...
    radia.simulink.solveHCurlEddyFosterHarmonic(model, 50, amplitude), amplitude);
verifyEqual(testCase, mean(lift(end - 2 * period + 1:end)), harmonic(3), "RelTol", 1e-2);
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
    plant + "/Moving HCurl Foster", "Parameters")), familyExpression);
verifyEqual(testCase, string(get_param( ...
    plant + "/Moving HCurl Foster", "FunctionName")), ...
    "radia_hcurl_eddy_foster_family_sfunction");
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
    "radia_maglev/MagLev Plant/Moving HCurl Foster", "FunctionName")), ...
    "radia_hcurl_eddy_foster_family_sfunction");
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
