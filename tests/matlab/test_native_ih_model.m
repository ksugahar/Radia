function tests = test_native_ih_model
%TEST_NATIVE_IH_MODEL Verify the tracked standalone native preview model.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
if bdIsLoaded("radia_simulink_library")
    close_system("radia_simulink_library",0);
end
clear radia_ih_eddy_sfun radia_ih_thermal_sfun
testCase.TestData.FileGenConfig = Simulink.fileGenControl("getConfig");
testCase.TestData.FileGenRoot = string(tempname("C:\temp"));
Simulink.fileGenControl("set", ...
    CacheFolder=fullfile(testCase.TestData.FileGenRoot,"cache"), ...
    CodeGenFolder=fullfile(testCase.TestData.FileGenRoot,"codegen"), ...
    createDir=true);
end

function teardownOnce(testCase)
Simulink.fileGenControl("setConfig", ...
    config=testCase.TestData.FileGenConfig);
if isfolder(testCase.TestData.FileGenRoot)
    rmdir(testCase.TestData.FileGenRoot,"s");
end
end

function testBuilderCreatesClosedNativeModel(testCase)
outputDirectory = string(tempname("C:\temp"));
mkdir(outputDirectory);
modelName = "radia_ih_model_" + erase(string(java.util.UUID.randomUUID),"-");
cleanup = onCleanup(@() closeIfLoaded(modelName));
modelPath = radia.simulink.buildIHNativeModel( ...
    ModelName=modelName, OutputDirectory=outputDirectory, Open=false);
verifyTrue(testCase,isfile(modelPath));
load_system(modelPath);
verifyEqual(testCase,string(get_param(modelName+"/Eddy","FunctionName")), ...
    "radia_ih_eddy_sfun");
verifyEqual(testCase,string(get_param(modelName+"/Thermal","FunctionName")), ...
    "radia_ih_thermal_sfun");
verifyEqual(testCase,string(get_param(modelName+"/IH Parameters","Mask")),"on");
contract=get_param(modelName+"/IH Parameters","UserData");
verifyEqual(testCase,string(contract.backend),"matlab-level2+radia-mex-handles");
verifyFalse(testCase,contract.python_fallback);
workspace=get_param(modelName,"ModelWorkspace");
config=workspace.getVariable("radia_ih_config");
verifyEqual(testCase,string(config.dt_order), ...
    "eddy;transport(theta_prev,theta_now);thermal");
verifyEqual(testCase,string(get_param(modelName,"Solver")),"FixedStepDiscrete");
set_param(modelName,"SimulationCommand","update");
sim(modelName,"StopTime","0.2");
end

function testTrackedModelLoadsAndUpdates(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
modelPath=fullfile(root,"matlab","radia_ih.slx");
verifyTrue(testCase,isfile(modelPath));
load_system(modelPath);
cleanup=onCleanup(@() closeIfLoaded("radia_ih"));
set_param("radia_ih","SimulationCommand","update");
verifyEqual(testCase,string(get_param("radia_ih/Eddy","FunctionName")), ...
    "radia_ih_eddy_sfun");
verifyEqual(testCase,string(get_param("radia_ih/Thermal","FunctionName")), ...
    "radia_ih_thermal_sfun");
end

function testPhysicalConfigRequiresVolReports(testCase)
modelName = "radia_ih_config_gate_" + ...
    erase(string(java.util.UUID.randomUUID),"-");
new_system(modelName);
cleanup = onCleanup(@() closeIfLoaded(modelName));
workspace = get_param(modelName,"ModelWorkspace");
config = radia.simulink.makeIHNativeSmokeConfig();
config = rmfield(config,"diagnostic_only");
config.operator_assembly = "preassembled";
workspace.assignin("radia_ih_config",config);
verifyError(testCase,@() radia.simulink.configureIHNativeModel(modelName), ...
    "radia:simulink:IHConfigVolCheck");
end

function testMakeConfigRequiresStrictLabelContract(testCase)
spec = nativeSpec();
verifyError(testCase,@() radia.simulink.makeIHNativeConfig( ...
    spec,NHeat=1,NTemperature=1,CellWeights=1,HeatCellWeights=1), ...
    "radia:simulink:IHConfigVolContract");
end

function testMakeConfigRejectsNonlinearPreview(testCase)
spec = nativeSpec();
spec.bh_mode = "nonlinear";
verifyError(testCase,@() radia.simulink.makeIHNativeConfig( ...
    spec,NHeat=1,NTemperature=1,CellWeights=1,HeatCellWeights=1), ...
    "radia:simulink:IHConfigBHMode");
end

function spec = nativeSpec()
spec = struct( ...
    "frequency",5e4,"current",100,"wp_vol","workpiece.vol", ...
    "peec_step","coil.step","method","peec-bem", ...
    "solver","direct","thermal_mesh_type","fem", ...
    "bh_mode","linear","n_eddy_unknown",1, ...
    "eddy_matrix_real",1,"eddy_matrix_imag",0, ...
    "eddy_rhs_real",1,"eddy_rhs_imag",0,"heat_projection",1, ...
    "mass_row_ptr",[0;1],"mass_col",0,"mass_value",1, ...
    "stiffness_row_ptr",[0;1],"stiffness_col",0, ...
    "stiffness_value",0,"initial_temperature_K",293.15);
end

function closeIfLoaded(name)
if bdIsLoaded(name), close_system(name,0); end
end
