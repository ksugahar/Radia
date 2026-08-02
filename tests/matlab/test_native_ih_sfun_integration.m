function tests = test_native_ih_sfun_integration
%TEST_NATIVE_IH_SFUN_INTEGRATION Exercise the native Eddy -> Thermal path.
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

function testHeatRaisesTemperature(testCase)
out1 = runNativeIH(1.0);
out2 = runNativeIH(2.0);
verifyGreaterThan(testCase,out1.temperature_K(end),293.15);
verifyEqual(testCase,out2.heat_density_W_per_m3(end), ...
    4*out1.heat_density_W_per_m3(end),"RelTol",1e-10);
verifyGreaterThan(testCase,out2.temperature_K(end),out1.temperature_K(end));
end

function out = runNativeIH(current)
model = "radia_native_ih_" + erase(string(java.util.UUID.randomUUID),"-");
cfg = radia.simulink.validateIHNativeConfig(contractConfig(1,1,1));
assignin("base","radia_native_ih_config",cfg);
new_system(model); cleanup=onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Constant",model+"/I","Value",num2str(current),"Position",[20 35 70 65]);
add_block("simulink/Sources/Constant",model+"/angle","Value","0","Position",[20 95 70 125]);
add_block("simulink/Sources/Constant",model+"/Tamb","Value","293.15","Position",[20 155 70 185]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Eddy", ...
    "FunctionName","radia_ih_eddy_sfun","Parameters","radia_native_ih_config", ...
    "Position",[120 45 240 105]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Thermal", ...
    "FunctionName","radia_ih_thermal_sfun","Parameters","radia_native_ih_config", ...
    "Position",[300 55 430 115]);
add_block("simulink/Sinks/To Workspace",model+"/HeatOut","VariableName","native_heat", ...
    "SaveFormat","Array","Position",[480 35 570 65]);
add_block("simulink/Sinks/To Workspace",model+"/TempOut","VariableName","native_temp", ...
    "SaveFormat","Array","Position",[480 100 570 130]);
add_line(model,"I/1","Eddy/1"); add_line(model,"angle/1","Eddy/2");
add_line(model,"Eddy/1","Thermal/1"); add_line(model,"Tamb/1","Thermal/2"); add_line(model,"angle/1","Thermal/3");
add_line(model,"Thermal/1","Eddy/3");
add_line(model,"Eddy/1","HeatOut/1"); add_line(model,"Thermal/1","TempOut/1");
set_param(model,"StopTime","0.3");
simOut=sim(model,"ReturnWorkspaceOutputs","on");
out.heat_density_W_per_m3=simOut.get("native_heat");
out.temperature_K=simOut.get("native_temp");
end

function testEddyAngleRotatesHeatInWorkpieceCoordinates(testCase)
heatAtOrigin = runEddyAtAngle(0);
heatAtQuarterTurn = runEddyAtAngle(pi/2);
verifyEqual(testCase,heatAtOrigin,[1,2,3,4],"AbsTol",1e-12);
verifyEqual(testCase,heatAtQuarterTurn,[2,3,4,1],"AbsTol",1e-12);
verifyEqual(testCase,sum(heatAtQuarterTurn),sum(heatAtOrigin), ...
    "AbsTol",1e-12);
end

function heat = runEddyAtAngle(angle)
model = "radia_native_ih_eddy_angle_" + ...
    erase(string(java.util.UUID.randomUUID),"-");
n = 4;
cfg = contractConfig(1,n,n);
cfg.heat_projection = (1:n).';
cfg = radia.simulink.validateIHNativeConfig(cfg);
assignin("base","radia_native_ih_eddy_angle_config",cfg);
new_system(model); cleanup=onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Constant",model+"/I","Value","1");
add_block("simulink/Sources/Constant",model+"/angle", ...
    "Value",compose("%.17g",angle));
add_block("simulink/Sources/Constant",model+"/T", ...
    "Value","293.15*ones(4,1)");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Eddy", ...
    "FunctionName","radia_ih_eddy_sfun", ...
    "Parameters","radia_native_ih_eddy_angle_config");
add_block("simulink/Sinks/To Workspace",model+"/HeatOut", ...
    "VariableName","rotated_heat","SaveFormat","Array");
add_line(model,"I/1","Eddy/1"); add_line(model,"angle/1","Eddy/2");
add_line(model,"T/1","Eddy/3"); add_line(model,"Eddy/1","HeatOut/1");
set_param(model,"StopTime","0","SolverType","Fixed-step", ...
    "Solver","FixedStepDiscrete","FixedStep","0.1");
simOut=sim(model,"ReturnWorkspaceOutputs","on");
values=simOut.get("rotated_heat");
heat=values(end,:);
end

function testRotationTransportConservesWeightedEnergy(testCase)
model = "radia_native_ih_rotation_" + erase(string(java.util.UUID.randomUUID),"-");
n = 4;
weights = [1;2;1;2];
initial = [300;340;380;420];
cfg = contractConfig(1,n,n);
cfg.heat_cell_weights = weights;
cfg.heat_to_temperature_projection = eye(n);
cfg.temperature_cell_weights = weights;
cfg.initial_temperature_K = initial;
cfg = radia.simulink.validateIHNativeConfig(cfg);
assignin("base","radia_native_ih_rotation_config",cfg);
new_system(model); cleanup=onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Constant",model+"/I","Value","0");
add_block("simulink/Sources/Constant",model+"/angle","Value","pi/2");
add_block("simulink/Sources/Constant",model+"/Tamb","Value","293.15");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Eddy", ...
    "FunctionName","radia_ih_eddy_sfun","Parameters", ...
    "radia_native_ih_rotation_config");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Thermal", ...
    "FunctionName","radia_ih_thermal_sfun","Parameters", ...
    "radia_native_ih_rotation_config");
add_block("simulink/Sinks/To Workspace",model+"/TempOut", ...
    "VariableName","rotation_temp","SaveFormat","Array");
add_line(model,"I/1","Eddy/1"); add_line(model,"angle/1","Eddy/2");
add_line(model,"Eddy/1","Thermal/1"); add_line(model,"Tamb/1","Thermal/2");
add_line(model,"angle/1","Thermal/3"); add_line(model,"Thermal/1","Eddy/3");
add_line(model,"Thermal/1","TempOut/1");
set_param(model,"StopTime","0.1","SolverType","Fixed-step", ...
    "Solver","FixedStepDiscrete","FixedStep","0.1");
simOut=sim(model,"ReturnWorkspaceOutputs","on");
temperature=simOut.get("rotation_temp");
verifyEqual(testCase,weights.'*temperature(end,:).',weights.'*initial, ...
    "AbsTol",1e-10);
verifyNotEqual(testCase,temperature(end,:),initial.');
end

function testThermalStateFeedsTemperatureDependentEddy(testCase)
cfg = contractConfig(1,1,1);
cfg.eddy_matrix_temperature_slope_real = 0.02;
cfg.eddy_matrix_temperature_slope_imag = 0;
cfg.bh_reference_temperature_K = 293.15;
cfg = radia.simulink.validateIHNativeConfig(cfg);
result = runClosedLoop(cfg,0.4);
heat = result.heat_density_W_per_m3(:,1);
temperature = result.temperature_K(:,1);
verifyGreaterThan(testCase,numel(heat),2);
verifyGreaterThan(testCase,heat(1),heat(end));
verifyGreaterThan(testCase,temperature(end),temperature(1));
verifyTrue(testCase,cfg.temperature_change_recomputes_eddy);
verifyFalse(testCase,cfg.current_change_recomputes_eddy);
end

function testConfigSupportsDistinctHeatAndTemperatureSpaces(testCase)
cfg = contractConfig(1,2,3);
cfg.heat_to_temperature_projection = [1,0;0.5,0.5;0,1];
cfg = radia.simulink.validateIHNativeConfig(cfg);
verifyEqual(testCase,numel(cfg.heat_to_temperature_projection),6);

bad = cfg;
bad.heat_to_temperature_projection = ones(2,2);
verifyError(testCase,@() radia.simulink.validateIHNativeConfig(bad), ...
    "radia:simulink:IHConfigMatrix");
end

function testConfigRejectsNonlinearAndInvalidCSR(testCase)
cfg = contractConfig(1,1,2);
cfg.bh_mode = "nonlinear";
verifyError(testCase,@() radia.simulink.validateIHNativeConfig(cfg), ...
    "radia:simulink:IHConfigBHMode");

cfg = contractConfig(1,1,2);
cfg.stiffness_col = flipud(cfg.stiffness_col);
verifyError(testCase,@() radia.simulink.validateIHNativeConfig(cfg), ...
    "radia:simulink:IHConfigThermalSparsity");
end

function testMatAndJsonConfigLoaders(testCase)
cfg = contractConfig(1,2,3);
matFile = string(tempname("C:\temp")) + ".mat";
jsonFile = string(tempname("C:\temp")) + ".json";
cleanup = onCleanup(@() deleteFiles(matFile,jsonFile));
config = cfg;
save(matFile,"config");
fid = fopen(jsonFile,"w");
fwrite(fid,jsonencode(cfg),"char");
fclose(fid);
fromMat = radia.simulink.validateIHNativeConfig( ...
    radia.simulink.loadIHNativeConfig(matFile));
fromJson = radia.simulink.validateIHNativeConfig( ...
    radia.simulink.loadIHNativeConfig(jsonFile));
verifyEqual(testCase,fromMat.n_temperature,3);
verifyEqual(testCase,fromJson.n_heat,2);
clear cleanup
deleteFiles(matFile,jsonFile);
end

function testSingularEddyFailsAndFreshLifecycleRecovers(testCase)
singular = contractConfig(1,1,1);
singular.eddy_matrix_real = 0;
singular = radia.simulink.validateIHNativeConfig(singular);
verifyError(testCase,@() runClosedLoop(singular,0), ...
    "Simulink:blocks:MSFB_BlockMethodFailed");

recovered = runClosedLoop( ...
    radia.simulink.validateIHNativeConfig(contractConfig(1,1,1)),0);
verifyTrue(testCase,all(isfinite(recovered.temperature_K),"all"));
end

function cfg = contractConfig(nUnknown,nHeat,nTemperature)
row = (0:nTemperature).';
column = (0:nTemperature-1).';
cfg = struct( ...
    "schema","radia.ih.simulink.native_sfunction.v1", ...
    "backend","matlab-level2+radia-mex-handles","python_fallback",false, ...
    "n_eddy_unknown",nUnknown,"n_heat",nHeat, ...
    "n_temperature",nTemperature,"bh_mode","linear", ...
    "eddy_solver","fem","thermal_solver","fem", ...
    "eddy_matrix_real",eye(nUnknown), ...
    "eddy_matrix_imag",zeros(nUnknown), ...
    "eddy_rhs_real",ones(nUnknown,1), ...
    "eddy_rhs_imag",zeros(nUnknown,1), ...
    "heat_projection",ones(nHeat,nUnknown), ...
    "heat_cell_weights",ones(nHeat,1), ...
    "heat_to_temperature_projection",ones(nTemperature,nHeat), ...
    "mass_row_ptr",row,"mass_col",column,"mass_value",ones(nTemperature,1), ...
    "stiffness_row_ptr",row,"stiffness_col",column, ...
    "stiffness_value",zeros(nTemperature,1), ...
    "temperature_cell_weights",ones(nTemperature,1), ...
    "initial_temperature_K",293.15*ones(nTemperature,1), ...
    "sample_time_s",0.1,"thermal_tolerance",1e-12, ...
    "thermal_max_iterations",100,"rotation_mode","periodic-uniform", ...
    "angle_origin_rad",0);
end

function result = runClosedLoop(cfg,stopTime)
model = "radia_native_ih_feedback_" + ...
    erase(string(java.util.UUID.randomUUID),"-");
assignin("base","radia_native_ih_feedback_config",cfg);
new_system(model); cleanup=onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Constant",model+"/I","Value","1");
add_block("simulink/Sources/Constant",model+"/angle","Value","0");
add_block("simulink/Sources/Constant",model+"/Tamb","Value","293.15");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Eddy", ...
    "FunctionName","radia_ih_eddy_sfun", ...
    "Parameters","radia_native_ih_feedback_config");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    model+"/Thermal", ...
    "FunctionName","radia_ih_thermal_sfun", ...
    "Parameters","radia_native_ih_feedback_config");
add_block("simulink/Sinks/To Workspace",model+"/HeatOut", ...
    "VariableName","feedback_heat","SaveFormat","Array");
add_block("simulink/Sinks/To Workspace",model+"/TempOut", ...
    "VariableName","feedback_temp","SaveFormat","Array");
add_line(model,"I/1","Eddy/1"); add_line(model,"angle/1","Eddy/2");
add_line(model,"Eddy/1","Thermal/1"); add_line(model,"Tamb/1","Thermal/2");
add_line(model,"angle/1","Thermal/3"); add_line(model,"Thermal/1","Eddy/3");
add_line(model,"Eddy/1","HeatOut/1"); add_line(model,"Thermal/1","TempOut/1");
set_param(model,"StopTime",num2str(stopTime,17), ...
    "SolverType","Fixed-step","Solver","FixedStepDiscrete", ...
    "FixedStep",num2str(cfg.sample_time_s,17));
simOut = sim(model,"ReturnWorkspaceOutputs","on");
result.heat_density_W_per_m3 = simOut.get("feedback_heat");
result.temperature_K = simOut.get("feedback_temp");
end

function deleteFiles(varargin)
for file = string(varargin)
    if isfile(file), delete(file); end
end
end

function closeIfLoaded(name)
if bdIsLoaded(name), close_system(name,0); end
end
