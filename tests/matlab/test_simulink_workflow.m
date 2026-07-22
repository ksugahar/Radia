function tests = test_simulink_workflow
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
end

function testIHPlantMatchesEnergyAndCoolingContract(testCase)
plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ...
    ThermalConductance_W_per_K=2, ...
    SampleTime_s=0.1, ...
    InitialTemperature_K=293.15);
time_s = (0:0.1:1).';
power_W = 100 * ones(size(time_s));
ambient_K = 293.15 * ones(size(time_s));
result = radia.simulink.simulateIHWaveform( ...
    plant, time_s, power_W, ambient_K);

verifyEqual(testCase, plant.schema, "radia.ih.simulink.plant.v1");
verifyEqual(testCase, result.schema, "radia.ih.simulink.waveform.v1");
verifyEqual(testCase, result.final_state(2), 110, "AbsTol", 1e-12);
verifyGreaterThan(testCase, result.temperature_K(end), 293.15);
verifyGreaterThanOrEqual(testCase, result.heat_loss_W(end), 0);
verifyEqual(testCase, result.position_rad, zeros(size(time_s)));
end

function testMotionAwarePowerLUTClipsOutsideTrainingRange(testCase)
lut = radia.simulink.makeIHPowerLUT( ...
    {[0, 1], [10, 20]}, [0, 10; 20, 30]);
power_W = radia.simulink.evaluateIHPowerLUT(lut, [0.5, 15; 2, 15]);
verifyEqual(testCase, power_W, [15; 25], "AbsTol", 1e-12);
verifyEqual(testCase, lut.schema, "radia.ih.simulink.power_lut.v1");
end

function testOneDimensionalPowerLUT(testCase)
lut = radia.simulink.makeIHPowerLUT({[0, 1]}, [5, 7]);
power_W = radia.simulink.evaluateIHPowerLUT(lut, [0.5; 2]);
verifyEqual(testCase, power_W, [6; 7], "AbsTol", 1e-12);
verifySize(testCase, lut.table_power_W, [2, 1]);
end

function testPeriodicEddyHeatDensityLUTTracksCurrentAndAngle(testCase)
lut = makeTestEddyLut();
heatDensity = radia.simulink.evaluateIHEddyHeatDensityLUT( ...
    lut, [10; 5; -5], [0; pi; 2 * pi + pi / 2]);
verifyEqual(testCase, lut.schema, ...
    "radia.ih.simulink.eddy_heat_density_lut.v1");
verifyEqual(testCase, heatDensity, [1000; 500; 375], "AbsTol", 1e-12);
verifyEqual(testCase, lut.rotation_angle_breakpoints_rad(end), 2 * pi, ...
    "AbsTol", 1e-12);
verifyEqual(testCase, lut.table_heat_density_W_per_m3(end, :, :), ...
    lut.table_heat_density_W_per_m3(1, :, :), "AbsTol", 1e-12);
end

function testDriveMotionAndTemperatureFeedback(testCase)
position = [0, 1];
drive = [0, 10];
temperature = [293.15, 303.15];
[p, d, t] = ndgrid(position, drive, temperature);
values = 1 + p + 0.1 * d + 0.01 * (t - 293.15);
lut = radia.simulink.makeIHPowerLUT( ...
    {position, drive, temperature}, values, ...
    InputNames=["position_rad", "drive_A", "temperature_K"]);
plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=0, SampleTime_s=0.1);
time_s = (0:0.1:0.3).';
result = radia.simulink.simulateIHDrive( ...
    plant, time_s, 5 * ones(size(time_s)), 293.15 * ones(size(time_s)), lut, ...
    Position_rad=0.5 * ones(size(time_s)), Speed_rad_s=ones(size(time_s)));
verifyEqual(testCase, result.power_in_W(1), 2, "AbsTol", 1e-12);
verifyGreaterThan(testCase, result.power_in_W(end), result.power_in_W(1));
verifyGreaterThan(testCase, result.final_state(1), 293.15);
verifyEqual(testCase, result.schema, "radia.ih.simulink.drive_waveform.v1");
end

function testOptimizationAdapter(testCase)
objective = @(x) (x(1) - 2)^2 + (x(2) + 1)^2;
result = radia.simulink.optimizeIH( ...
    objective, [0, 0], [-3, -3], [3, 3], ...
    UseFmincon=false, Display="off", MaxIterations=200);
verifyEqual(testCase, result.schema, "radia.ih.simulink.optimization.v1");
verifyLessThan(testCase, norm(result.x - [2, -1]), 1e-3);
verifyLessThan(testCase, result.objective_value, 1e-6);
end

function testSimulinkBuilderUpdatesPlantAndPIDModels(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
eddyLut = makeTestEddyLut();
for includePID = [false, true]
    modelName = "radia_ih_test_model_" + string(includePID);
    cleanup = onCleanup(@() closeIfLoaded(modelName));
    radia.simulink.buildIHControlModel(modelName, plant, eddyLut, ...
        IncludePID=includePID, Save=false, Open=false);
    set_param(modelName, "SimulationCommand", "update");
    verifyTrue(testCase, bdIsLoaded(modelName));
    eddyPath = modelName + "/Eddy Current";
    thermalPath = modelName + "/Thermal";
    parameterPath = modelName + "/IH Parameters";
    verifyEqual(testCase, string(get_param(eddyPath, "BlockType")), "SubSystem");
    verifyEqual(testCase, string(get_param(thermalPath, "BlockType")), "SubSystem");
    verifyEqual(testCase, string(get_param(parameterPath, "BlockType")), "SubSystem");
    verifyEqual(testCase, string(get_param(eddyPath, "Mask")), "on");
    verifyEqual(testCase, string(get_param(thermalPath, "Mask")), "on");
    verifyEqual(testCase, string(get_param(parameterPath, "Mask")), "on");
    parameterMask = Simulink.Mask.get(parameterPath);
    verifyEqual(testCase, ...
        str2double(parameterMask.getParameter("angle_period_rad").Value), ...
        eddyLut.angle_period_rad, "AbsTol", eps(eddyLut.angle_period_rad));
    verifyEqual(testCase, ...
        eval(parameterMask.getParameter("region_volumes_m3").Value), ...
        eddyLut.region_volumes_m3);
    verifyEqual(testCase, ...
        string(get_param(eddyPath + "/HeatDensity_1", "BlockType")), ...
        "Lookup_n-D");
    verifyEqual(testCase, ...
        string(get_param(eddyPath + "/PreviousRotationAngle", "BlockType")), ...
        "UnitDelay");
    currentMagnitudePorts = get_param(eddyPath + "/CurrentMagnitude", "PortHandles");
    verifyNotEqual(testCase, get_param(currentMagnitudePorts.Inport, "Line"), -1);
    verifyEqual(testCase, ...
        string(get_param(thermalPath + "/ThermalPlant", "BlockType")), ...
        "DiscreteStateSpace");
    eddyPorts = get_param(eddyPath, "PortHandles");
    thermalPorts = get_param(thermalPath, "PortHandles");
    parameterPorts = get_param(parameterPath, "PortHandles");
    eddyLine = get_param(eddyPorts.Outport(1), "Line");
    verifyTrue(testCase, any(get_param(eddyLine, "DstPortHandle") == ...
        thermalPorts.Inport(1)));
    verifyTrue(testCase, any(get_param( ...
        get_param(parameterPorts.Outport(1), "Line"), "DstPortHandle") == ...
        eddyPorts.Inport(3)));
    verifyTrue(testCase, any(get_param( ...
        get_param(parameterPorts.Outport(2), "Line"), "DstPortHandle") == ...
        eddyPorts.Inport(4)));
    verifyTrue(testCase, any(get_param( ...
        get_param(parameterPorts.Outport(3), "Line"), "DstPortHandle") == ...
        thermalPorts.Inport(3)));
    if includePID
        delayPath = modelName + "/TemperatureFeedbackDelay";
        verifyEqual(testCase, string(get_param(delayPath, "BlockType")), ...
            "UnitDelay");
        verifyEqual(testCase, str2double(get_param(delayPath, "InitialCondition")), ...
            plant.x0(1), "AbsTol", eps(plant.x0(1)));
    end
    clear cleanup
    closeIfLoaded(modelName);
end
end

function testOpenIHLaunchEntryPoint(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

target = radia.simulink.openIH(Open=false);
verifyEqual(testCase, target, ...
    "radia_simulink_library/Applications/Induction Heating");
verifyTrue(testCase, bdIsLoaded("radia_simulink_library"));
verifyEqual(testCase, string(get_param(target, "Mask")), "on");
close_system("radia_simulink_library", 0);
end

function testSavedIHObjectUsesRequestedOutputDirectory(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
eddyLut = makeTestEddyLut();
modelName = "radia_ih_saved_object_test";
outputDirectory = string(fullfile("C:\temp", modelName));
modelPath = fullfile(outputDirectory, modelName + ".slx");
cleanup = onCleanup(@() cleanSavedModel(modelName, outputDirectory));
actualPath = radia.simulink.openIH( ...
    ModelName=modelName, Plant=plant, EddyLUT=eddyLut, ...
    OutputDirectory=outputDirectory, Save=true, Open=false);

verifyEqual(testCase, actualPath, modelPath);
verifyTrue(testCase, isfile(modelPath));
verifyEqual(testCase, string(get_param(modelName + "/Eddy Current", "Mask")), "on");
verifyEqual(testCase, string(get_param(modelName + "/Thermal", "Mask")), "on");
verifyEqual(testCase, string(get_param(modelName + "/IH Parameters", "Mask")), "on");
clear cleanup
cleanSavedModel(modelName, outputDirectory);
end

function testIHSampleObjectHasSourcesAndSeparatedSubsystems(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

modelName = "radia_ih_sample_object_test";
outputDirectory = string(fullfile("C:\temp", modelName));
cleanup = onCleanup(@() cleanSavedModel(modelName, outputDirectory));
modelPath = radia.simulink.buildIHSampleModel( ...
    ModelName=modelName, OutputDirectory=outputDirectory, Open=false);

verifyTrue(testCase, isfile(modelPath));
verifyEqual(testCase, string(get_param(modelName + "/coil_current_rms_A", ...
    "BlockType")), "Step");
verifyEqual(testCase, string(get_param(modelName + "/rotation_angle_rad", ...
    "ReferenceBlock")), "simulink/Sources/Ramp");
verifyEqual(testCase, string(get_param(modelName + "/ambient_temperature_K", ...
    "BlockType")), "Constant");
verifyEqual(testCase, string(get_param(modelName + "/Eddy Current", "Mask")), "on");
verifyEqual(testCase, string(get_param(modelName + "/Thermal", "Mask")), "on");
verifyEqual(testCase, string(get_param(modelName + "/IH Parameters", "Mask")), "on");
verifyEqual(testCase, string(get_param(modelName + "/Thermal/ThermalPlant", ...
    "BlockType")), "DiscreteStateSpace");
verifyEmpty(testCase, find_system(modelName, "FollowLinks", "on", ...
    "LookUnderMasks", "all", "BlockType", "S-Function"));
set_param(modelName, "SimulationCommand", "update");
clear cleanup
cleanSavedModel(modelName, outputDirectory);
end

function testTrackedIHSampleArtifactIsMatlabOnly(testCase)
hasSimulink = exist("load_system", "file") == 2 || ...
    exist("load_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
sampleFile = fullfile(repoRoot, "matlab", "radia_ih_sample.slx");
verifyTrue(testCase, isfile(sampleFile));
cleanup = onCleanup(@() closeIfLoaded("radia_ih_sample"));
load_system(sampleFile);
set_param("radia_ih_sample", "SimulationCommand", "update");

verifyEqual(testCase, string(get_param( ...
    "radia_ih_sample/rotation_angle_rad", "ReferenceBlock")), ...
    "simulink/Sources/Ramp");
verifyEqual(testCase, string(get_param( ...
    "radia_ih_sample/Eddy Current", "Mask")), "on");
verifyEqual(testCase, string(get_param( ...
    "radia_ih_sample/Thermal", "Mask")), "on");
verifyEmpty(testCase, find_system("radia_ih_sample", "FollowLinks", "on", ...
    "LookUnderMasks", "all", "BlockType", "S-Function"));
clear cleanup
closeIfLoaded("radia_ih_sample");
end

function testEddyAngleUsesPreviousSample(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=0, SampleTime_s=0.1);
eddyLut = makeTestEddyLut();
modelName = "radia_ih_angle_delay_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildIHControlModel( ...
    modelName, plant, eddyLut, StopTime_s=0.3, Save=false, Open=false);
time_s = (0:0.1:0.3).';
angle = [0; pi / 2; pi; 3 * pi / 2];
inputData = [time_s, 10 * ones(size(time_s)), angle, ...
    293.15 * ones(size(time_s))];
assignin("base", "radia_ih_angle_delay_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_ih_angle_delay_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
heatDensity = dataset.getElement(5).Values.Data;
verifyEqual(testCase, heatDensity, [1000; 1000; 1500; 2000], ...
    "AbsTol", 1e-10);
clear cleanup
closeIfLoaded(modelName);
end

function testSimulinkExternalInputProducesHeatingWaveform(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1, ...
    InitialTemperature_K=293.15);
modelName = "radia_ih_external_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
eddyLut = makeTestEddyLut();
radia.simulink.buildIHControlModel( ...
    modelName, plant, eddyLut, Save=false, Open=false);
time_s = (0:0.1:1).';
inputData = [time_s, 10 * ones(size(time_s)), zeros(size(time_s)), ...
    293.15 * ones(size(time_s))];
assignin("base", "radia_ih_external_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_ih_external_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
temperatureSignal = dataset.getElement(1);
heatDensitySignal = dataset.getElement(5);
verifyGreaterThan(testCase, temperatureSignal.Values.Data(end), 293.15);
verifyEqual(testCase, heatDensitySignal.Values.Data(end), 1000, "AbsTol", 1e-12);
clear cleanup
closeIfLoaded(modelName);
end

function testRadiaSFunctionBlockMatchesStandardPlant(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
modelName = "radia_ih_sfunction_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
eddyLut = makeTestEddyLut();
radia.simulink.buildIHControlModel(modelName, plant, eddyLut, ...
    PlantBlock="radia-sfunction", StopTime_s=1.0, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.1:1).';
inputData = [time_s, 10 * ones(size(time_s)), zeros(size(time_s)), ...
    293.15 * ones(size(time_s))];
assignin("base", "radia_ih_sfunction_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_ih_sfunction_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
temperatureSignal = dataset.getElement(1);
reference = radia.simulink.simulateIHWaveform( ...
    plant, time_s, 100 * ones(size(time_s)), 293.15 * ones(size(time_s)));
verifyGreaterThan(testCase, temperatureSignal.Values.Data(end), 293.15);
verifyEqual(testCase, temperatureSignal.Values.Data(end), ...
    reference.temperature_K(end), "AbsTol", 1e-10);
clear cleanup
closeIfLoaded(modelName);
end

function testTeam28CLNLUTAndSimulinkBlock(testCase)
lut = radia.simulink.makeTeam28CLNLUT();
[force_N, lift_N, slope_N_per_m] = radia.simulink.evaluateTeam28CLNForce( ...
    lut, 0.0, 20.0);
verifyEqual(testCase, lut.schema, "radia.team28.cln_lut.v1");
verifyEqual(testCase, lut.frequency_Hz, 50.0);
verifyEqual(testCase, force_N, -1.096266057556509, "AbsTol", 1e-12);
verifyEqual(testCase, lift_N, -force_N, "AbsTol", 1e-12);
verifyGreaterThan(testCase, slope_N_per_m, 0);

quarterForce = radia.simulink.evaluateTeam28CLNForce(lut, 0.0, 10.0);
verifyEqual(testCase, quarterForce, force_N / 4, "AbsTol", 1e-12);

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    return
end
modelName = "radia_team28_cln_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildTeam28CLNModel(modelName, lut, ...
    SampleTime_s=0.1, StopTime_s=0.2, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.1:0.2).';
inputData = [time_s, zeros(size(time_s)), 20 * ones(size(time_s))];
assignin("base", "radia_team28_cln_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_team28_cln_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
forceSignal = dataset.getElement(1);
verifyEqual(testCase, forceSignal.Values.Data(end), force_N, "AbsTol", 1e-12);
clear cleanup
closeIfLoaded(modelName);
end

function testNativeMexStateSpaceSimulinkBlocks(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
modelName = "radia_ih_native_mex_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
eddyLut = makeTestEddyLut();
radia.simulink.buildIHControlModel(modelName, plant, eddyLut, ...
    PlantBlock="radia-mex", StopTime_s=1.0, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.1:1).';
inputData = [time_s, 10 * ones(size(time_s)), zeros(size(time_s)), ...
    293.15 * ones(size(time_s))];
assignin("base", "radia_ih_native_mex_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_ih_native_mex_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
temperatureSignal = dataset.getElement(1);
reference = radia.simulink.simulateIHWaveform( ...
    plant, time_s, 100 * ones(size(time_s)), 293.15 * ones(size(time_s)));
verifyEqual(testCase, temperatureSignal.Values.Data(end), ...
    reference.temperature_K(end), "AbsTol", 1e-10);
clear cleanup
closeIfLoaded(modelName);

R = [2.0, 0.1; 0.1, 1.0];
L = [3.0, 0.2; 0.2, 2.0];
P = [1.0; 0.5];
stateModel = radia.simulink.makeHCurlEddyCLNModel( ...
    R, L, P, SampleTime_s=0.01);
modelName = "radia_hcurl_native_mex_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildHCurlEddyCLNModel(modelName, stateModel, ...
    Block="radia-mex", StopTime_s=0.02, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.01:0.02).';
assignin("base", "radia_hcurl_native_mex_input", [time_s, ones(size(time_s))]);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_hcurl_native_mex_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
responseSignal = dataset.getElement(1);
expectedState = stateModel.Bd + stateModel.Ad * stateModel.Bd;
verifyEqual(testCase, responseSignal.Values.Data(end), ...
    stateModel.Cd * expectedState, "AbsTol", 1e-12);
clear cleanup
closeIfLoaded(modelName);
end

function testHCurlEddyCLNStateSpaceAndMexSolve(testCase)
R = [2.0, 0.1; 0.1, 1.0];
L = [3.0, 0.2; 0.2, 2.0];
P = [1.0; 0.5];
model = radia.simulink.makeHCurlEddyCLNModel( ...
    R, L, P, SampleTime_s=0.01);
verifyEqual(testCase, model.schema, "radia.hcurl.eddy_cln.state_space.v1");
verifyTrue(testCase, model.passive);
verifyEqual(testCase, size(model.A), [2, 2]);
verifyEqual(testCase, size(model.C), [1, 2]);

frequency_Hz = 50.0;
coefficients = radia.simulink.solveHCurlEddyCLNHarmonic(model, frequency_Hz, 2.0);
s = 1i * 2.0 * pi * frequency_Hz;
reference = (R + s * L) \ (-s * P * 2.0);
verifyEqual(testCase, coefficients, reference, "AbsTol", 1e-12);

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    return
end
modelName = "radia_hcurl_eddy_cln_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildHCurlEddyCLNModel(modelName, model, ...
    StopTime_s=0.02, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.01:0.02).';
inputData = [time_s, ones(size(time_s))];
assignin("base", "radia_hcurl_eddy_cln_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_hcurl_eddy_cln_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
responseSignal = dataset.getElement(1);
expectedState = model.Bd + model.Ad * model.Bd;
verifyEqual(testCase, responseSignal.Values.Data(end), ...
    model.Cd * expectedState, "AbsTol", 1e-12);
clear cleanup
closeIfLoaded(modelName);
end

function testHCurlEddyCLNExchangeLoadAndForce(testCase)
payload = struct( ...
    "schema", "radia.hcurl.eddy_cln.exchange.v1", ...
    "state_order", 2, ...
    "port_count", 1, ...
    "sample_time_s", 0.01, ...
    "has_sibc_termination", false, ...
    "arrays", struct( ...
        "resistance", struct("shape", [2, 2], "values", [2, 0.1, 0.1, 1]), ...
        "inductance", struct("shape", [2, 2], "values", [3, 0.2, 0.2, 2]), ...
        "surface_mass", struct("shape", [2, 2], "values", [0, 0, 0, 0]), ...
        "port_rhs", struct("shape", [2, 1], "values", [1, 0.5]), ...
        "force_operator", struct("shape", [3, 2, 1], ...
            "values", [1, 2, 3, 4, 5, 6])), ...
    "metadata", struct("frequency_hz", 50));
fileName = fullfile(tempdir, "radia_hcurl_exchange_test.json");
cleanup = onCleanup(@() deleteIfExists(fileName));
fid = fopen(fileName, "w");
fwrite(fid, jsonencode(payload), "char");
fclose(fid);

model = radia.simulink.loadHCurlEddyCLNModel(fileName);
verifyEqual(testCase, model.exchange_schema, ...
    "radia.hcurl.eddy_cln.exchange.v1");
verifyEqual(testCase, model.resistance, [2, 0.1; 0.1, 1], "AbsTol", 0);
verifyEqual(testCase, model.metadata.frequency_hz, 50);
force_N = radia.simulink.evaluateHCurlEddyCLNForce(model, [2; 3], 4);
verifyEqual(testCase, force_N, [16; 36; 56], "AbsTol", 1e-12);
clear cleanup
deleteIfExists(fileName);
end

function testHCurlEddyCLNHeightFamilyAndMovingBlock(testCase)
snapshot0 = makeFamilySnapshot(-1.0, 1.0, [1; 2; 3]);
snapshot1 = makeFamilySnapshot(1.0, 3.0, [3; 4; 5]);
payload = struct( ...
    "schema", "radia.hcurl.eddy_cln.family.v1", ...
    "shared_state_basis", true, ...
    "sample_time_s", 0.01, ...
    "state_order", 1, ...
    "port_count", 1, ...
    "snapshots", [snapshot0, snapshot1], ...
    "metadata", struct("frequency_hz", 50));
fileName = fullfile(tempdir, "radia_hcurl_family_test.json");
cleanupFile = onCleanup(@() deleteIfExists(fileName));
fid = fopen(fileName, "w");
fwrite(fid, jsonencode(payload), "char");
fclose(fid);

family = radia.simulink.loadHCurlEddyCLNFamily(fileName);
verifyEqual(testCase, family.snapshot_count, 2);
verifyEqual(testCase, family.positions_m, [-1; 1], "AbsTol", 0);
mid = radia.simulink.interpolateHCurlEddyCLNFamily(family, 0.0, ...
    BuildStateSpace=true);
verifyEqual(testCase, mid.resistance, 2.0, "AbsTol", 1e-12);
verifyEqual(testCase, mid.force_operator, [2; 3; 4], "AbsTol", 1e-12);
force_N = radia.simulink.evaluateHCurlEddyCLNForce(mid, 2.0, 4.0);
verifyEqual(testCase, force_N, [8; 12; 16], "AbsTol", 1e-12);
verifyError(testCase, ...
    @() radia.simulink.interpolateHCurlEddyCLNFamily(family, 2.0), ...
    "radia:simulink:HCurlCLNExtrapolation");

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if hasSimulink
    modelName = "radia_hcurl_eddy_cln_family_test";
    cleanupModel = onCleanup(@() closeIfLoaded(modelName));
    radia.simulink.buildHCurlEddyCLNFamilyModel(modelName, family, ...
        StopTime_s=0.02, Save=false, Open=false);
    set_param(modelName, "SimulationCommand", "update");
    verifyTrue(testCase, bdIsLoaded(modelName));
    clear cleanupModel
    closeIfLoaded(modelName);
end
clear cleanupFile
deleteIfExists(fileName);
end

function testRadiaLibraryBrowserArtifact(testCase)
output="C:\temp\radia_simulink_library_test";
library=radia.simulink.buildLibrary(OutputDirectory=output);
verifyTrue(testCase,isfile(library)); load_system(library); cleanup=onCleanup(@()closeIfLoaded("radia_simulink_library"));
verifyEqual(testCase,string(get_param("radia_simulink_library","BlockDiagramType")),"library");
applications = [ ...
    "Electromagnet", "em"; ...
    "PCB PEEC", "pcb"; ...
    "Motor", "motor"; ...
    "Stream Function", "streamfunction"; ...
    "Induction Heating", "ih"];
for row = 1:size(applications, 1)
    blockPath = "radia_simulink_library/Applications/" + applications(row, 1);
    verifyEqual(testCase, string(get_param(blockPath, "FunctionName")), ...
        "radia_application_sfun");
    verifyEqual(testCase, string(get_param(blockPath, "Mask")), "on");
    parameters = string(get_param(blockPath, "Parameters"));
    verifyTrue(testCase, startsWith(parameters, "'" + applications(row, 2) + "'"));
    displayScript = string(get_param(blockPath, "MaskDisplay"));
    verifyTrue(testCase, contains(displayScript, "'run'"));
    verifyTrue(testCase, contains(displayScript, "'status'"));
    verifyTrue(testCase, contains(displayScript, "'primary'"));
    verifyTrue(testCase, contains(displayScript, "'elapsed_s'"));
end
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/LTspice Circuit","FunctionName")),"radia_ltspice_sfun");
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/Hysteretic LTspice Plant","FunctionName")),"radia_hysteretic_ltspice_sfun");
verifyEqual(testCase,string(get_param("radia_simulink_library/Optimization/Optuna Optimization","FunctionName")),"radia_optuna_sfun");
clear cleanup
end

function testApplicationBlockCompilesInModel(testCase)
output = "C:\temp\radia_simulink_application_compile_test";
library = radia.simulink.buildLibrary(OutputDirectory=output);
load_system(library);
libraryCleanup = onCleanup(@() closeIfLoaded("radia_simulink_library"));

modelName = "radia_application_compile_test";
closeIfLoaded(modelName);
new_system(modelName);
modelCleanup = onCleanup(@() closeIfLoaded(modelName));
add_block("simulink/Sources/Constant", modelName + "/Trigger", ...
    Value="false", OutDataTypeStr="boolean");
add_block("radia_simulink_library/Applications/Electromagnet", ...
    modelName + "/Application");
add_line(modelName, "Trigger/1", "Application/1");
set_param(modelName, "SimulationCommand", "update");

verifyTrue(testCase, bdIsLoaded(modelName));
verifyEqual(testCase, ...
    string(get_param(modelName + "/Application", "FunctionName")), ...
    "radia_application_sfun");
clear modelCleanup libraryCleanup
end

function testApplicationConfigAndFailureArtifacts(testCase)
root = "C:\temp\radia_simulink_application_test";
config = fullfile(root, "em_config.json");
radia.simulink.writeApplicationConfig("em", struct(), config, ...
    PrimaryKey="B_origin_mag_T");
decoded = jsondecode(fileread(config));
verifyEqual(testCase, string(decoded.schema), ...
    "radia.simulink.application_config.v1");
verifyEqual(testCase, string(decoded.application), "em");

result = radia.simulink.runApplication("em", config, ...
    RunRoot=root, Timeout_s=10, ThrowOnFailure=false);
verifyEqual(testCase, string(result.status), "failed");
verifyTrue(testCase, contains(string(result.error), "Coil script"));
verifyTrue(testCase, isfile(result.result_json));
verifyTrue(testCase, isfile(result.log));
end

function testApplicationLauncherFailureKeepsArtifacts(testCase)
root = "C:\temp\radia_simulink_launcher_failure_test";
config = fullfile(root, "em_config.json");
radia.simulink.writeApplicationConfig("em", struct(), config);

result = radia.simulink.runApplication("em", config, ...
    RunRoot=root, ...
    PythonExecutable=fullfile(root, "missing_python.exe"), ...
    ThrowOnFailure=false);

verifyEqual(testCase, string(result.status), "failed");
verifyNotEqual(testCase, result.process_status, 0);
verifyTrue(testCase, contains(string(result.error), ...
    "did not create result.json"));
verifyTrue(testCase, isfile(result.result_json));
verifyTrue(testCase, isfile(result.log));
verifyTrue(testCase, isfile(result.command_file));
decoded = jsondecode(fileread(result.result_json));
verifyEqual(testCase, string(decoded.radia_result.schema), ...
    "radia.simulink.application_run.v1");
verifyEqual(testCase, string(decoded.radia_result.status), "failed");
end

function lut = makeTestEddyLut()
theta = [0; pi / 2; pi; 3 * pi / 2];
current = [0; 5; 10];
angleScale = [1; 1.5; 2; 1.5];
currentScale = (current / 10).^2;
heatDensity = 1000 * angleScale * currentScale.';
lut = radia.simulink.makeIHEddyHeatDensityLUT( ...
    theta, current, heatDensity, RegionVolumes_m3=0.1, ...
    CarrierFrequency_Hz=50e3, Source="test map");
end

function snapshot = makeFamilySnapshot(height_m, resistance, force_operator)
base = radia.simulink.makeHCurlEddyCLNModel( ...
    resistance, 1.0, 1.0, SampleTime_s=0.01);
snapshot = struct( ...
    "height_m", height_m, ...
    "arrays", struct( ...
        "resistance", familyArray([1, 1], resistance), ...
        "inductance", familyArray([1, 1], 1.0), ...
        "surface_mass", familyArray([1, 1], 0.0), ...
        "port_rhs", familyArray([1, 1], 1.0), ...
        "force_operator", familyArray([3, 1, 1], force_operator)), ...
    "metadata", struct("height_offset_m", height_m), ...
    "state_order", 1, ...
    "port_count", 1, ...
    "sample_time_s", 0.01);
if isempty(base)
    error("radia:test:unreachable", "base model construction failed.");
end
end

function encoded = familyArray(shape, values)
encoded = struct("shape", shape, "values", values(:).');
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end

function deleteIfExists(fileName)
if isfile(fileName)
    delete(fileName);
end
end

function cleanSavedModel(modelName, outputDirectory)
closeIfLoaded(modelName);
if isfolder(outputDirectory)
    rmdir(outputDirectory, "s");
end
end
