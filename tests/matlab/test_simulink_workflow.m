function tests = test_simulink_workflow
tests = functiontests(localfunctions);
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
for includePID = [false, true]
    modelName = "radia_ih_test_model_" + string(includePID);
    cleanup = onCleanup(@() closeIfLoaded(modelName));
    radia.simulink.buildIHControlModel(modelName, plant, ...
        IncludePID=includePID, Save=false, Open=false);
    set_param(modelName, "SimulationCommand", "update");
    verifyTrue(testCase, bdIsLoaded(modelName));
    clear cleanup
    closeIfLoaded(modelName);
end
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
radia.simulink.buildIHControlModel(modelName, plant, Save=false, Open=false);
time_s = (0:0.1:1).';
inputData = [time_s, 100 * ones(size(time_s)), 293.15 * ones(size(time_s))];
assignin("base", "radia_ih_external_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_ih_external_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simOut = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simOut.get("yout");
temperatureSignal = dataset.getElement(1);
verifyGreaterThan(testCase, temperatureSignal.Values.Data(end), 293.15);
clear cleanup
closeIfLoaded(modelName);
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
