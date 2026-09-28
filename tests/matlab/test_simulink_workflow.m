function tests = test_simulink_workflow
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.RepoRoot = string(repoRoot);
scratch = string(getenv("RADIA_TEMP_ROOT"));
if strlength(scratch) == 0
    if ispc, scratch = "C:\temp"; else, scratch = string(tempdir); end
end
if ~isfolder(scratch), mkdir(scratch); end
root = string(tempname(scratch));
mkdir(root);
testCase.TestData.ScratchRoot = root;
testCase.addTeardown(@() rmdir(root, "s"));
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

function testAdjointTopologyOptimizationBlock(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

runner = radia.topopt.makeAdjointDemoRunner();
modelName = "radia_adjoint_topopt_test";
cleanup = onCleanup(@() cleanupAdjointModel(modelName));
blockPath = radia.simulink.buildAdjointOptimizationBlock(modelName, ...
    Runner=runner,Solver="mma",Save=false);
add_block("simulink/Sources/Constant",modelName + "/Start", ...
    Value="1");
signalNames = ["Objective","Status","Iterations","Design"];
variables = ["objective_signal","status_signal", ...
    "iterations_signal","design_signal"];
for index = 1:numel(signalNames)
    add_block("simulink/Sinks/To Workspace", ...
        modelName + "/" + signalNames(index), ...
        VariableName=variables(index),SaveFormat="Array");
    add_line(modelName, ...
        "Adjoint Topology Optimization/" + index, ...
        signalNames(index) + "/1");
end
add_line(modelName,"Start/1","Adjoint Topology Optimization/1");
set_param(modelName,"StopTime","0.2","Solver","FixedStepDiscrete", ...
    "FixedStep","0.1","ReturnWorkspaceOutputs","on");
simulation = sim(modelName);

objectiveSignal = simulation.get("objective_signal");
statusSignal = simulation.get("status_signal");
iterationsSignal = simulation.get("iterations_signal");
designSignal = simulation.get("design_signal");
verifyEqual(testCase,size(objectiveSignal,2),1);
verifyEqual(testCase,size(statusSignal,2),1);
verifyEqual(testCase,size(iterationsSignal,2),1);
verifyEqual(testCase,size(designSignal,2),2);
if hasOptimizationToolbox()
    verifyEqual(testCase,statusSignal(end),2);
    verifyTrue(testCase,runner.Result.converged,runner.Result.output.message);
    verifyLessThan(testCase,runner.Result.objective,0.006);
    verifyLessThanOrEqual(testCase,max(runner.Result.constraints),1e-6);
    verifyEqual(testCase,objectiveSignal(end),runner.Result.objective, ...
        "AbsTol",1e-12);
    verifyEqual(testCase,designSignal(end,:).',runner.Result.design, ...
        "AbsTol",1e-12);
    verifyEqual(testCase,iterationsSignal(end),height(runner.Result.history)-1);
else
    verifyEqual(testCase,statusSignal(end),-1);
    verifyTrue(testCase,isnan(objectiveSignal(end)));
    verifyEqual(testCase,iterationsSignal(end),0);
    verifyEqual(testCase,designSignal(end,:).',runner.InitialDesign, ...
        "AbsTol",0);
end
verifyEqual(testCase,string(get_param(blockPath,"FunctionName")), ...
    "radia_adjoint_optimization_sfun");
clear cleanup
cleanupAdjointModel(modelName);
end

function testAdjointLibraryBlockCompilesInModel(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

output = fullfile(testCase.TestData.ScratchRoot, "adjoint_library");
library = radia.simulink.buildLibrary(OutputDirectory=output);
load_system(library);
libraryCleanup = onCleanup(@() closeIfLoaded("radia_simulink_library"));
modelName = "radia_adjoint_library_compile_test";
modelCleanup = onCleanup(@() cleanupAdjointModel(modelName));
new_system(modelName);
add_block("simulink/Sources/Constant",modelName + "/Start",Value="0");
add_block( ...
    "radia_simulink_library/Optimization/Adjoint Topology Optimization", ...
    modelName + "/Topology Optimization");
add_line(modelName,"Start/1","Topology Optimization/1");
set_param(modelName,"SimulationCommand","update");

blockPath = modelName + "/Topology Optimization";
verifyEqual(testCase,string(get_param(blockPath,"FunctionName")), ...
    "radia_adjoint_optimization_sfun");
verifyEqual(testCase,string(get_param(blockPath,"Parameters")), ...
    "runner,sample_time_s");
verifyTrue(testCase,isa(evalin("base","radia_adjoint_runner"), ...
    "radia.topopt.AdjointRunner"));
clear modelCleanup libraryCleanup
cleanupAdjointModel(modelName);
closeIfLoaded("radia_simulink_library");
end

function testStreamFunctionOptimizationModel(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

modelName = "radia_streamfunction_topopt_test";
output = fullfile(testCase.TestData.ScratchRoot, "streamfunction");
modelPath = fullfile(output,modelName + ".slx");
cleanup = onCleanup(@() cleanupStreamFunctionModel(modelName,modelPath));
radia.simulink.buildStreamFunctionOptimizationModel( ...
    ModelName=modelName,OutputDirectory=output, ...
    SampleTime_s=0.1,HistoryCapacity=64);
load_system(modelPath);
set_param(modelName + "/Run Adjoint","Value","true");
ports = [4,5,6,7,8,9,10,11];
variables = [ ...
    "sf_objective","sf_status","sf_iterations","sf_design", ...
    "sf_history_count","sf_history_iteration", ...
    "sf_history_objective","sf_history_constraint"];
for index = 1:numel(ports)
    add_block("simulink/Sinks/To Workspace", ...
        modelName + "/" + variables(index), ...
        VariableName=variables(index),SaveFormat="Array");
    add_line(modelName, ...
        "Stream Function Optimization/" + ports(index), ...
        variables(index) + "/1");
end
set_param(modelName,"StopTime","0.2","ReturnWorkspaceOutputs","on");
simulation = sim(modelName);

objective = simulation.get("sf_objective");
status = simulation.get("sf_status");
iterations = simulation.get("sf_iterations");
design = simulation.get("sf_design");
historyCount = simulation.get("sf_history_count");
historyIteration = simulation.get("sf_history_iteration");
historyObjective = simulation.get("sf_history_objective");
historyConstraint = simulation.get("sf_history_constraint");
count = historyCount(end);
verifyEqual(testCase,size(design,2),2);
verifyEqual(testCase,size(historyIteration,2),64);
verifyEqual(testCase,size(historyObjective,2),64);
verifyEqual(testCase,size(historyConstraint,2),64);
if hasOptimizationToolbox()
    verifyEqual(testCase,status(end),2);
    verifyLessThan(testCase,objective(end),0.003);
    verifyGreaterThan(testCase,iterations(end),1);
    verifyGreaterThan(testCase,count,2);
    verifyTrue(testCase,all(isfinite(historyObjective(end,1:count))));
    verifyLessThanOrEqual(testCase,historyConstraint(end,count),1e-6);
else
    verifyEqual(testCase,status(end),-1);
    verifyTrue(testCase,isnan(objective(end)));
    verifyEqual(testCase,iterations(end),0);
    verifyEqual(testCase,count,0);
end
blockPath = modelName + "/Stream Function Optimization";
contract = get_param(blockPath,"UserData");
verifyFalse(testCase,contract.python_per_step);
verifyEqual(testCase,string(contract.adjoint_backend), ...
    "matlab-native-analytic-gradient");
customRunner = radia.topopt.makeStreamFunctionAdjointDemoRunner();
assignin("base","radia_streamfunction_custom_runner",customRunner);
set_param(blockPath,"radia_streamfunction_adjoint_runner", ...
    "radia_streamfunction_custom_runner");
set_param(modelName,"SimulationCommand","update");
verifyEqual(testCase,string(get_param( ...
    blockPath + "/Analytic Topology Optimization","Parameters")), ...
    "radia_streamfunction_custom_runner,0.10000000000000001,64");
clear cleanup
cleanupStreamFunctionModel(modelName,modelPath);
end

function available = hasOptimizationToolbox()
available = ~isempty(ver("optim")) && ...
    license("test","Optimization_Toolbox") && ...
    exist("fmincon","file") == 2;
end

function testTrackedStreamFunctionOptimizationArtifact(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

testDirectory = fileparts(mfilename("fullpath"));
repositoryRoot = fileparts(fileparts(testDirectory));
modelPath = fullfile(repositoryRoot,"matlab", ...
    "radia_streamfunction_optimization.slx");
verifyTrue(testCase,isfile(modelPath));
modelName = "radia_streamfunction_optimization";
closeIfLoaded(modelName);
cleanup = onCleanup(@() closeIfLoaded(modelName));
load_system(modelPath);
set_param(modelName,"SimulationCommand","update");
blockPath = modelName + "/Stream Function Optimization";
verifyEqual(testCase,string(get_param(blockPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    blockPath + "/Analytic Topology Optimization","FunctionName")), ...
    "radia_streamfunction_topology_sfun");
contract = get_param(blockPath,"UserData");
verifyEqual(testCase,string(contract.domain), ...
    "stream-function-topology-optimization");
verifyFalse(testCase,contract.python_per_step);
clear cleanup
closeIfLoaded(modelName);
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

function testTeam28LiftLUTAndSimulinkBlock(testCase)
lut = radia.simulink.makeTeam28LiftLUT();
[force_N, lift_N, slope_N_per_m] = radia.simulink.evaluateTeam28LiftForce( ...
    lut, 0.0, 20.0);
verifyEqual(testCase, lut.schema, "radia.team28.lift_lut.v1");
verifyEqual(testCase, lut.frequency_Hz, 50.0);
verifyTrue(testCase, lut.source_validation_passed);
sweep = jsondecode(fileread(fullfile(testCase.TestData.RepoRoot, ...
    "validation_test", "maglev", "team28_axisym_fem_height_sweep.json")));
verifyEqual(testCase, lift_N, sweep.upward_physical_N(sweep.dZ_mm == 0), ...
    "AbsTol", 1e-12);
verifyEqual(testCase, force_N, -lift_N, "AbsTol", 1e-12);
verifyGreaterThan(testCase, slope_N_per_m, 0);

quarterForce = radia.simulink.evaluateTeam28LiftForce(lut, 0.0, 10.0);
verifyEqual(testCase, quarterForce, force_N / 4, "AbsTol", 1e-12);

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    return
end
modelName = "radia_team28_lift_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildTeam28LiftModel(modelName, lut, ...
    SampleTime_s=0.1, StopTime_s=0.2, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.1:0.2).';
inputData = [time_s, zeros(size(time_s)), 20 * ones(size(time_s))];
assignin("base", "radia_team28_lift_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_team28_lift_input", ...
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

stateModel = radia.simulink.makeHCurlEddyFosterModel( ...
    [0.5; 2.0], [1.0; 0.5], SampleTime_s=0.01);
modelName = "radia_hcurl_native_mex_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildHCurlEddyFosterModel(modelName, stateModel, ...
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

function testHCurlEddyFosterStateSpaceAndHarmonicSolve(testCase)
% The Foster modes of an independent R/L/P model reproduce its dense solve.
R = [2.0, 0.1; 0.1, 1.0];
L = [3.0, 0.2; 0.2, 2.0];
P = [1.0; 0.5];
[V, Lambda] = eig(R, L);
V = V ./ sqrt(diag(V.' * L * V)).';
model = radia.simulink.makeHCurlEddyFosterModel( ...
    diag(Lambda), V.' * P, SampleTime_s=0.01);
verifyEqual(testCase, model.schema, "radia.hcurl.eddy_foster.state_space.v1");
verifyTrue(testCase, model.passive);
verifyEqual(testCase, size(model.A), [2, 2]);
verifyEqual(testCase, size(model.C), [1, 2]);
continuous = expm([model.A, model.B; zeros(1, 3)] * model.sample_time_s);
verifyEqual(testCase, model.Ad, continuous(1:2, 1:2), "AbsTol", 1e-14);
verifyEqual(testCase, model.Bd, continuous(1:2, 3), "AbsTol", 1e-14);

frequency_Hz = 50.0;
z = radia.simulink.solveHCurlEddyFosterHarmonic(model, frequency_Hz, 2.0);
s = 1i * 2.0 * pi * frequency_Hz;
reference = (R + s * L) \ (-s * P * 2.0);
verifyEqual(testCase, V * z, reference, "AbsTol", 1e-12);

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    return
end
modelName = "radia_hcurl_eddy_foster_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildHCurlEddyFosterModel(modelName, model, ...
    StopTime_s=0.02, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:0.01:0.02).';
inputData = [time_s, ones(size(time_s))];
assignin("base", "radia_hcurl_eddy_foster_input", inputData);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_hcurl_eddy_foster_input", ...
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

function testHCurlEddyFosterExchangeLoadAndForce(testCase)
payload = struct( ...
    "schema", "radia.hcurl.eddy_foster.exchange.v1", ...
    "state_order", 2, ...
    "port_count", 1, ...
    "sample_time_s", 0.01, ...
    "has_sibc_termination", false, ...
    "decay_rates", familyArray(2, [2.0, 7.0]), ...
    "arrays", struct( ...
        "modal_port_rhs", familyArray([2, 1], [1, 0.5]), ...
        "modal_force_operator", familyArray([3, 2, 1], [1, 2, 3, 4, 5, 6])), ...
    "metadata", struct("frequency_hz", 50));
fileName = fullfile(tempdir, "radia_hcurl_foster_exchange_test.json");
cleanup = onCleanup(@() deleteIfExists(fileName));
fid = fopen(fileName, "w");
fwrite(fid, jsonencode(payload), "char");
fclose(fid);

model = radia.simulink.loadHCurlEddyFosterModel(fileName);
verifyEqual(testCase, model.exchange_schema, ...
    "radia.hcurl.eddy_foster.exchange.v1");
verifyEqual(testCase, model.decay_rates, [2.0; 7.0], "AbsTol", 0);
verifyEqual(testCase, model.modal_port_rhs, [1; 0.5], "AbsTol", 0);
verifyEqual(testCase, model.metadata.frequency_hz, 50);
force_N = radia.simulink.evaluateHCurlEddyFosterForce(model, [2; 3], 4);
verifyEqual(testCase, force_N, [16; 36; 56], "AbsTol", 1e-12);
clear cleanup
deleteIfExists(fileName);
end

function testHCurlEddyFosterHeightFamilyAndMovingBlock(testCase)
payload = struct( ...
    "schema", "radia.hcurl.eddy_foster.family.v1", ...
    "shared_modes", true, ...
    "sample_time_s", 0.01, ...
    "state_order", 1, ...
    "port_count", 1, ...
    "decay_rates", familyArray(1, 2.0), ...
    "snapshots", [makeFamilySnapshot(-1.0, 1.0, [1; 2; 3]), ...
        makeFamilySnapshot(1.0, 3.0, [3; 4; 5])], ...
    "metadata", struct("frequency_hz", 50));
fileName = fullfile(tempdir, "radia_hcurl_foster_family_test.json");
cleanupFile = onCleanup(@() deleteIfExists(fileName));
fid = fopen(fileName, "w");
fwrite(fid, jsonencode(payload), "char");
fclose(fid);

family = radia.simulink.loadHCurlEddyFosterFamily(fileName);
verifyEqual(testCase, family.snapshot_count, 2);
verifyEqual(testCase, family.positions_m, [-1; 1], "AbsTol", 0);
mid = radia.simulink.interpolateHCurlEddyFosterFamily(family, 0.0);
verifyEqual(testCase, mid.modal_port_rhs, 2.0, "AbsTol", 1e-12);
verifyEqual(testCase, mid.Ad, exp(-2.0 * 0.01), "AbsTol", 1e-15);
verifyEqual(testCase, mid.force_operator, [2; 3; 4], "AbsTol", 1e-12);
force_N = radia.simulink.evaluateHCurlEddyFosterForce(mid, 2.0, 4.0);
verifyEqual(testCase, force_N, [8; 12; 16], "AbsTol", 1e-12);
verifyError(testCase, ...
    @() radia.simulink.interpolateHCurlEddyFosterFamily(family, 2.0), ...
    "radia:simulink:HCurlFosterExtrapolation");

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if hasSimulink
    modelName = "radia_hcurl_eddy_foster_family_test";
    cleanupModel = onCleanup(@() closeIfLoaded(modelName));
    radia.simulink.buildHCurlEddyFosterFamilyModel(modelName, family, ...
        StopTime_s=0.02, Save=false, Open=false);
    set_param(modelName, "SimulationCommand", "update");
    verifyTrue(testCase, bdIsLoaded(modelName));
    clear cleanupModel
    closeIfLoaded(modelName);
end
clear cleanupFile
deleteIfExists(fileName);
end

function testNativePeriodicMotorAngleFamilyBlock(testCase)
grid = [0; pi];
A = reshape([1.0, 0.5], 1, 1, 2);
B = reshape([1.0, 3.0], 1, 1, 2);
C = reshape([1.0, 3.0], 1, 1, 2);
D = zeros(1, 1, 2);
Q = reshape([2.0, 4.0], 1, 1, 2);
R = reshape([1.0, 2.0], 1, 1, 2);
S = reshape([0.0, 2.0], 1, 1, 2);
family = radia.simulink.makeMotorAngleFamily( ...
    grid, A, B, C, D, Q, R, S, 2.0, ...
    Period_rad=2*pi, SampleTime_s=0.01);
verifyEqual(testCase, family.schema, "radia.motor.periodic-angle-family.v1");
verifyEqual(testCase, family.backend, "native-mex-periodic-interpolation");
verifyEqual(testCase, family.output_count, 2);

hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    return
end
modelName = "radia_motor_angle_family_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildMotorAngleFamilyModel(modelName, family, ...
    StopTime_s=0.02, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
block = modelName + "/MotorAngleFamily";
verifyEqual(testCase, string(get_param(block, "FunctionName")), ...
    "radia_motor_angle_family_mex_sfunction");
workspace = get_param(modelName, "ModelWorkspace");
contract = getVariable(workspace, "radia_motor_angle_family_contract");
verifyFalse(testCase, contract.python_per_step);
verifyFalse(testCase, contract.matlab_matrix_algebra_per_step);

time = (0:family.sample_time_s:0.02).';
inputs = Simulink.SimulationData.Dataset;
inputs = inputs.addElement(timeseries( ...
    repmat(pi/2, numel(time), 1), time), "mechanical_angle_rad");
inputs = inputs.addElement(timeseries( ...
    repmat(3.0, numel(time), 1), time), "model_inputs");
simulation = Simulink.SimulationInput(char(modelName));
simulation = simulation.setExternalInput(inputs);
simulation = simulation.setModelParameter( ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simulationOutput = sim(simulation);
logged = simulationOutput.get("yout").getElement(1).Values.Data;
expected = [ ...
    4.0, 19.5; ...
    15.0, 122.625; ...
    23.25, 259.5234375];
verifyEqual(testCase, logged, expected, "AbsTol", 1e-10);
clear cleanup
closeIfLoaded(modelName);
end

function testRadiaLibraryBrowserArtifact(testCase)
output=fullfile(testCase.TestData.ScratchRoot, "library");
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
ihMask=Simulink.Mask.get(ihPath);
ihConfigParameter=ihMask.getParameter("config_file");
verifyEqual(testCase,string(ihConfigParameter.Evaluate),"off");
verifyTrue(testCase,contains(string(ihConfigParameter.Callback), ...
    "onIHConfigFileChanged"));
verifyEqual(testCase,string(get_param( ...
    ihPath+"/Eddy Level-2 S-Function","FunctionName")), ...
    "radia_ih_eddy_sfun");
verifyEqual(testCase,string(get_param( ...
    ihPath+"/Thermal Level-2 S-Function","FunctionName")), ...
    "radia_ih_thermal_sfun");
ihPorts=get_param(ihPath,"PortHandles");
verifyEqual(testCase,numel(ihPorts.Inport),3);
verifyEqual(testCase,numel(ihPorts.Outport),3);
verifyEqual(testCase,string(get_param( ...
    ihPath+"/IH Monitor/Bus","OutDataTypeStr")),"Bus: IHMonitorBusV1");
verifyEqual(testCase,string(get_param( ...
    ihPath+"/IH Monitor/Runtime","FunctionName")), ...
    "radia_ih_monitor_sfun");
verifyEqual(testCase,string(get_param(ihPath+"/IH Monitor","Mask")),"on");
verifyEmpty(testCase,find_system(ihPath,SearchDepth=1, ...
    BlockType="BusCreator"));
verifyEmpty(testCase,find_system(ihPath, ...
    "LookUnderMasks","all","BlockType","Lookup_n-D"));
verifyEmpty(testCase,find_system(ihPath, ...
    "LookUnderMasks","all","BlockType","DiscreteStateSpace"));
ihContract=get_param(ihPath,"UserData");
verifyEqual(testCase,string(ihContract.backend),"matlab-level2+radia-mex-handles");
verifyFalse(testCase,ihContract.python_fallback);
verifyTrue(testCase,ihContract.distributed_field);
verifyFalse(testCase,ihContract.surrogate);
verifyEqual(testCase,string(ihContract.monitor_bus),"IHMonitorBusV1");
verifyEqual(testCase,string(ihContract.temperature_feedback), ...
    "previous-accepted-thermal-state");
maglevPath="radia_simulink_library/Applications/Magnetic Levitation";
verifyEqual(testCase,string(get_param(maglevPath,"BlockType")),"SubSystem");
verifyEqual(testCase,string(get_param(maglevPath,"Mask")),"on");
maglevMask=Simulink.Mask.get(maglevPath);
maglevFamilyParameter=maglevMask.getParameter("family");
verifyEqual(testCase,string(maglevFamilyParameter.Evaluate),"off");
verifyTrue(testCase,contains(string(maglevFamilyParameter.Callback), ...
    "onMagLevBlockFamilyChanged"));
maglevPorts=get_param(maglevPath,"PortHandles");
verifyEqual(testCase,numel(maglevPorts.Inport),3);
verifyEqual(testCase,numel(maglevPorts.Outport),2);
verifyEqual(testCase,string(get_param( ...
    maglevPath+"/Moving HCurl Foster","FunctionName")), ...
    "radia_hcurl_eddy_foster_family_sfunction");
maglevContract=get_param(maglevPath,"UserData");
verifyEqual(testCase,string(maglevContract.backend), ...
    "matlab-level2-shared-mode-foster");
verifyFalse(testCase,maglevContract.python_per_step);
verifyFalse(testCase,maglevContract.surrogate);
fieldStatsPath="radia_simulink_library/Utilities/Field Stats";
verifyEqual(testCase,string(get_param(fieldStatsPath,"BlockType")), ...
    "SubSystem");
verifyEqual(testCase,string(get_param(fieldStatsPath,"Tag")), ...
    "RadiaFieldStats");
verifyEqual(testCase,string(get_param(fieldStatsPath,"Mask")),"on");
fieldStatsPorts=get_param(fieldStatsPath,"PortHandles");
verifyEqual(testCase,numel(fieldStatsPorts.Inport),1);
verifyEqual(testCase,numel(fieldStatsPorts.Outport),1);
sfOptimizationPath = ...
    "radia_simulink_library/Applications/Stream Function Optimization";
verifyEqual(testCase,string(get_param(sfOptimizationPath,"BlockType")), ...
    "SubSystem");
verifyEqual(testCase,string(get_param(sfOptimizationPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    sfOptimizationPath + "/Stream Function Design","FunctionName")), ...
    "radia_application_sfun");
verifyEqual(testCase,string(get_param( ...
    sfOptimizationPath + "/Analytic Topology Optimization", ...
    "FunctionName")),"radia_streamfunction_topology_sfun");
sfPorts = get_param(sfOptimizationPath,"PortHandles");
verifyEqual(testCase,numel(sfPorts.Inport),2);
verifyEqual(testCase,numel(sfPorts.Outport),11);
sfContract = get_param(sfOptimizationPath,"UserData");
verifyEqual(testCase,string(sfContract.domain), ...
    "stream-function-topology-optimization");
verifyFalse(testCase,sfContract.python_per_step);
verifyEqual(testCase,string(sfContract.visualization),"simulink-scope-xy");
verifyGreaterThan(testCase,getSimulinkBlockHandle( ...
    sfOptimizationPath + "/Objective History"),0);
verifyGreaterThan(testCase,getSimulinkBlockHandle( ...
    sfOptimizationPath + "/Constraint History"),0);
sfOptunaPath = ...
    "radia_simulink_library/Applications/Stream Function Optuna";
verifyEqual(testCase,string(get_param(sfOptunaPath,"FunctionName")), ...
    "radia_streamfunction_optuna_sfun");
verifyEqual(testCase,string(get_param(sfOptunaPath,"Mask")),"on");
sfOptunaPorts = get_param(sfOptunaPath,"PortHandles");
verifyEqual(testCase,numel(sfOptunaPorts.Inport),1);
verifyEqual(testCase,numel(sfOptunaPorts.Outport),7);
sfOptunaContract = get_param(sfOptunaPath,"UserData");
verifyEqual(testCase,string(sfOptunaContract.optimizer_backend), ...
    "matlab-native-optuna");
verifyEqual(testCase,string(sfOptunaContract.inner_factorization), ...
    "aca_plus_qr_tsvd");
verifyTrue(testCase,sfOptunaContract.python_per_trial);
verifyFalse(testCase,sfOptunaContract.python_per_step);
emOptimizationPath = ...
    "radia_simulink_library/Applications/Electromagnet Topology Optimization";
verifyEqual(testCase,string(get_param(emOptimizationPath,"BlockType")), ...
    "SubSystem");
verifyEqual(testCase,string(get_param(emOptimizationPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    emOptimizationPath + "/Electromagnet Analysis","FunctionName")), ...
    "radia_application_sfun");
verifyEqual(testCase,string(get_param( ...
    emOptimizationPath + "/Density Topology Optimization", ...
    "FunctionName")),"radia_electromagnet_topology_sfun");
emPorts = get_param(emOptimizationPath,"PortHandles");
verifyEqual(testCase,numel(emPorts.Inport),2);
verifyEqual(testCase,numel(emPorts.Outport),12);
emContract = get_param(emOptimizationPath,"UserData");
verifyEqual(testCase,string(emContract.domain), ...
    "electromagnet-topology-optimization");
verifyEqual(testCase,emContract.state_solves_per_gradient,1);
verifyEqual(testCase,emContract.adjoint_solves_per_gradient,1);
verifyFalse(testCase,emContract.finite_difference_fallback);
verifyFalse(testCase,emContract.python_per_optimization_step);
motorFamilyPath = ...
    "radia_simulink_library/Reduced Models/Motor Angle Family";
verifyEqual(testCase,string(get_param(motorFamilyPath,"FunctionName")), ...
    "radia_motor_angle_family_mex_sfunction");
verifyEqual(testCase,string(get_param(motorFamilyPath,"Mask")),"on");
verifyEqual(testCase,string(get_param(motorFamilyPath,"Parameters")), ...
    "family");
motorFamilyContract=get_param(motorFamilyPath,"UserData");
verifyEqual(testCase,string(motorFamilyContract.backend), ...
    "native-mex-periodic-interpolation");
verifyEqual(testCase,string(motorFamilyContract.state_lifecycle), ...
    "outputs-read;update-advance");
verifyFalse(testCase,motorFamilyContract.python_per_step);
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/LTspice Circuit","FunctionName")),"radia_ltspice_sfun");
verifyEqual(testCase,string(get_param("radia_simulink_library/LTspice/Hysteretic LTspice Plant","FunctionName")),"radia_hysteretic_ltspice_sfun");
optunaPath = "radia_simulink_library/Optimization/Optuna Study";
verifyEqual(testCase,string(get_param(optunaPath,"BlockType")),"SubSystem");
verifyEqual(testCase,string(get_param(optunaPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    optunaPath, ...
    "sampler_name")),"auto");
optunaPorts = get_param(optunaPath,"PortHandles");
verifyEqual(testCase,numel(optunaPorts.Inport),2);
verifyEqual(testCase,numel(optunaPorts.Outport),5);
verifyEqual(testCase,string(get_param( ...
    optunaPath + "/Monitor Bus","OutDataTypeStr")), ...
    "Bus: OptunaMonitorBusV1");
verifyEqual(testCase,string(get_param( ...
    optunaPath + "/Advanced Runtime","FunctionName")), ...
    "radia_optuna_sfun");
adjointPath = ...
    "radia_simulink_library/Optimization/Adjoint Topology Optimization";
verifyEqual(testCase,string(get_param(adjointPath,"FunctionName")), ...
    "radia_adjoint_optimization_sfun");
verifyEqual(testCase,string(get_param(adjointPath,"Mask")),"on");
adjointContract = get_param(adjointPath,"UserData");
verifyEqual(testCase,string(adjointContract.domain), ...
    "topology-optimization");
verifyEqual(testCase,string(adjointContract.backend), ...
    "matlab-native-adjoint");
verifyFalse(testCase,adjointContract.finite_difference_fallback);
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
output = fullfile(testCase.TestData.ScratchRoot, "application_compile");
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
add_block("radia_simulink_library/Applications/Stream Function", ...
    modelName + "/Stream Function");
add_line(modelName, "Trigger/1", "Application/1");
add_line(modelName, "Trigger/1", "Stream Function/1");
set_param(modelName, "SimulationCommand", "update");

verifyTrue(testCase, bdIsLoaded(modelName));
verifyEqual(testCase, ...
    string(get_param(modelName + "/Application", "FunctionName")), ...
    "radia_application_sfun");
verifyEqual(testCase, ...
    string(get_param(modelName + "/Stream Function", "FunctionName")), ...
    "radia_application_sfun");
clear modelCleanup libraryCleanup
end

function testApplicationConfigAndFailureArtifacts(testCase)
root = fullfile(testCase.TestData.ScratchRoot, "application");
config = fullfile(root, "em_config.json");
radia.simulink.writeApplicationConfig("em", struct(), config, ...
    PrimaryKey="B_origin_mag_T");
decoded = jsondecode(fileread(config));
verifyEqual(testCase, string(decoded.schema), ...
    "radia.simulink.application_config.v1");
verifyEqual(testCase, string(decoded.application), "em");

streamConfig = fullfile(root, "streamfunction_config.json");
streamSettings = struct( ...
    "method", "Design", ...
    "coil_vol", "coil.vol", ...
    "eval_vol", "dsv.vol", ...
    "target_cf", "1", ...
    "aca_eps", 2.5e-9);
radia.simulink.writeApplicationConfig( ...
    "streamfunction", streamSettings, streamConfig);
streamDecoded = jsondecode(fileread(streamConfig));
verifyEqual(testCase, string(streamDecoded.application), "streamfunction");
verifyEqual(testCase, streamDecoded.settings.aca_eps, 2.5e-9, ...
    "AbsTol", eps(2.5e-9));

result = radia.simulink.runApplication("em", config, ...
    RunRoot=root, Timeout_s=10, ThrowOnFailure=false);
verifyEqual(testCase, string(result.status), "failed");
verifyTrue(testCase, contains(string(result.error), "Coil script"));
verifyTrue(testCase, isfile(result.result_json));
verifyTrue(testCase, isfile(result.log));
end

function testApplicationLauncherFailureKeepsArtifacts(testCase)
root = fullfile(testCase.TestData.ScratchRoot, "launcher_failure");
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

function snapshot = makeFamilySnapshot(height_m, modal_port, force_operator)
snapshot = struct( ...
    "height_m", height_m, ...
    "arrays", struct( ...
        "modal_port_rhs", familyArray([1, 1], modal_port), ...
        "modal_force_operator", familyArray([3, 1, 1], force_operator)), ...
    "metadata", struct("height_offset_m", height_m), ...
    "state_order", 1, ...
    "port_count", 1, ...
    "sample_time_s", 0.01);
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

function cleanupAdjointModel(modelName)
closeIfLoaded(modelName);
evalin("base","clear radia_adjoint_runner");
end

function cleanupStreamFunctionModel(modelName,modelPath)
closeIfLoaded(modelName);
if isfile(modelPath)
    delete(modelPath);
end
evalin('base', ...
    'clear radia_streamfunction_adjoint_runner radia_streamfunction_custom_runner');
end
