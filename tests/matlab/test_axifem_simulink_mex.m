function tests = test_axifem_simulink_mex
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"), "-begin");
radia.setup(Force=true);
end

function testQ2StateModelContract(testCase)
mu0 = 4*pi*1e-7;
model = radia.simulink.makeAxiEddyElementModel( ...
    0, 1e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7, SampleTime_s=1e-5);
verifyEqual(testCase, string(model.schema), ...
    "radia.axifem.q2_eddy.state_space.v1");
verifyEqual(testCase, string(model.backend), "shared-native-q2-mex");
verifyTrue(testCase, model.element.axis_touching);
verifyEqual(testCase, model.free_dofs, 9);
verifyEqual(testCase, model.state_order, 1);
verifyFalse(testCase, model.python_per_step);
verifyLessThan(testCase, model.A, 0);
verifyGreaterThan(testCase, model.Ad, 0);
verifyLessThan(testCase, model.Ad, 1);
verifyEqual(testCase, model.Cd(1:8), zeros(8, 1), AbsTol=0);
verifyEqual(testCase, model.Cd(9), 1, AbsTol=0);
end

function testAxisDofsCannotBeReleased(testCase)
mu0 = 4*pi*1e-7;
call = @() radia.simulink.makeAxiEddyElementModel( ...
    0, 1e-3, 0, 1e-3, mu0, 5.8e7, DirichletDofs=[1, 2, 3]);
verifyError(testCase, call, "radia:simulink:AxiEddyAxis");
end

function testNativeMexSimulinkRun(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

mu0 = 4*pi*1e-7;
sampleTime = 1e-5;
model = radia.simulink.makeAxiEddyElementModel( ...
    0, 1e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7, SampleTime_s=sampleTime);
modelName = "radia_axifem_q2_eddy_mex_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
radia.simulink.buildAxiEddyElementModel(modelName, model, ...
    StopTime_s=2*sampleTime, Save=false, Open=false);
set_param(modelName, "SimulationCommand", "update");
time_s = (0:sampleTime:2*sampleTime).';
assignin("base", "radia_axifem_q2_eddy_input", [time_s, ones(size(time_s))]);
set_param(modelName, "LoadExternalInput", "on", ...
    "ExternalInput", "radia_axifem_q2_eddy_input", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
simulation = sim(modelName, "ReturnWorkspaceOutputs", "on");
dataset = simulation.get("yout");
signal = dataset.getElement(1);
expectedState = model.Bd + model.Ad * model.Bd;
verifyEqual(testCase, signal.Values.Data(end, :).', ...
    model.Cd * expectedState, "AbsTol", 1e-12);
modelWorkspace = get_param(modelName, "ModelWorkspace");
contract = evalin(modelWorkspace, "radia_axifem_eddy_contract");
verifyEqual(testCase, string(contract.backend), "native-mex-sfunction");
verifyFalse(testCase, contract.python_per_step);
clear cleanup
closeIfLoaded(modelName);
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
