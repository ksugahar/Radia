function tests = test_simulink_workflow
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
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

function testOpenIHLaunchEntryPoint(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

target = radia.simulink.openIH(Open=false);
verifyEqual(testCase, target, "radia_ih");
verifyTrue(testCase, bdIsLoaded("radia_ih"));
verifyEqual(testCase, string(get_param(target + "/Eddy", ...
    "FunctionName")), "radia_ih_eddy_sfun");
verifyEqual(testCase, string(get_param(target + "/Thermal", ...
    "FunctionName")), "radia_ih_thermal_sfun");
verifyEqual(testCase, string(get_param(target + "/IH Parameters", ...
    "Mask")), "on");
close_system("radia_ih", 0);
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

function testNativeMexHCurlStateSpaceSimulinkBlock(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

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
    "Stream Function", "streamfunction"];
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
ihPath="radia_simulink_library/Applications/Induction Heating";
verifyEqual(testCase,string(get_param(ihPath,"BlockType")),"SubSystem");
verifyEqual(testCase,string(get_param(ihPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    ihPath+"/Eddy MEX S-Function","FunctionName")), ...
    "radia_ih_eddy_sfun");
verifyEqual(testCase,string(get_param( ...
    ihPath+"/Thermal MEX S-Function","FunctionName")), ...
    "radia_ih_thermal_sfun");
ihPorts=get_param(ihPath,"PortHandles");
verifyEqual(testCase,numel(ihPorts.Inport),3);
verifyEqual(testCase,numel(ihPorts.Outport),2);
verifyEmpty(testCase,find_system(ihPath, ...
    "LookUnderMasks","all","BlockType","Lookup_n-D"));
verifyEmpty(testCase,find_system(ihPath, ...
    "LookUnderMasks","all","BlockType","DiscreteStateSpace"));
ihContract=get_param(ihPath,"UserData");
verifyEqual(testCase,string(ihContract.backend),"native-mex-sfunction");
verifyFalse(testCase,ihContract.python_fallback);
verifyTrue(testCase,ihContract.distributed_field);
verifyFalse(testCase,ihContract.surrogate);
verifyEqual(testCase,string(ihContract.temperature_feedback), ...
    "previous-accepted-thermal-state");
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/LTspice Circuit","FunctionName")),"radia_ltspice_sfun");
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/Hysteretic LTspice Plant","FunctionName")),"radia_hysteretic_ltspice_sfun");
verifyEqual(testCase,string(get_param("radia_simulink_library/Optimization/Optuna Optimization","FunctionName")),"radia_optuna_sfun");
verifyEqual(testCase,string(get_param( ...
    "radia_simulink_library/Optimization/Optuna Optimization","Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    "radia_simulink_library/Optimization/Optuna Optimization", ...
    "sampler_name")),"auto");
sheetPath="radia_simulink_library/Optimization/Sheet Metal Optimization";
verifyEqual(testCase,string(get_param(sheetPath,"FunctionName")),"radia_optuna_sfun");
verifyEqual(testCase,string(get_param(sheetPath,"Mask")),"on");
verifyEqual(testCase,string(get_param(sheetPath,"sampler_name")),"auto");
sheetContract=get_param(sheetPath,"UserData");
verifyEqual(testCase,string(sheetContract.domain),"sheet-metal");
verifyEqual(testCase,string(sheetContract.backend), ...
    "matlab-native-ngsolve-cubit");
verifyFalse(testCase,sheetContract.browser_required);
monitorPath="radia_simulink_library/Optimization/Optuna Monitor";
verifyEqual(testCase,string(get_param(monitorPath,"Mask")),"on");
monitorContract=get_param(monitorPath,"UserData");
verifyFalse(testCase,monitorContract.browser_required);
verifyEqual(testCase,string(monitorContract.visualization),"simulink-scope-xy");
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
