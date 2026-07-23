function tests = test_native_ih_sfun_integration
%TEST_NATIVE_IH_SFUN_INTEGRATION Exercise the native Eddy -> Thermal path.
tests = functiontests(localfunctions);
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
cfg = struct("n_eddy_unknown",1,"n_heat",1,"n_temperature",1, ...
    "bh_mode",'linear', ...
    "eddy_matrix_real",1,"eddy_matrix_imag",0, ...
    "eddy_rhs_real",1,"eddy_rhs_imag",0,"heat_projection",1, ...
    "mass_row_ptr",[0;1],"mass_col",0,"mass_value",1, ...
    "stiffness_row_ptr",[0;1],"stiffness_col",0,"stiffness_value",0, ...
    "temperature_cell_weights",1,"initial_temperature_K",293.15, ...
    "sample_time_s",0.1,"thermal_tolerance",1e-12, ...
    "thermal_max_iterations",100);
assignin("base","radia_native_ih_config",cfg);
new_system(model); cleanup=onCleanup(@() closeIfLoaded(model)); %#ok<NASGU>
add_block("simulink/Sources/Constant",model+"/I","Value",num2str(current),"Position",[20 35 70 65]);
add_block("simulink/Sources/Constant",model+"/angle","Value","0","Position",[20 95 70 125]);
add_block("simulink/Sources/Constant",model+"/T","Value","293.15","Position",[20 155 70 185]);
add_block("simulink/User-Defined Functions/S-Function",model+"/Eddy", ...
    "FunctionName","radia_ih_eddy_sfun","Parameters","radia_native_ih_config", ...
    "Position",[120 45 240 105]);
add_block("simulink/User-Defined Functions/S-Function",model+"/Thermal", ...
    "FunctionName","radia_ih_thermal_sfun","Parameters","radia_native_ih_config", ...
    "Position",[300 55 430 115]);
add_block("simulink/Sinks/To Workspace",model+"/HeatOut","VariableName","native_heat", ...
    "SaveFormat","Array","Position",[480 35 570 65]);
add_block("simulink/Sinks/To Workspace",model+"/TempOut","VariableName","native_temp", ...
    "SaveFormat","Array","Position",[480 100 570 130]);
add_line(model,"I/1","Eddy/1"); add_line(model,"angle/1","Eddy/2"); add_line(model,"T/1","Eddy/3");
add_line(model,"Eddy/1","Thermal/1"); add_line(model,"T/1","Thermal/2"); add_line(model,"angle/1","Thermal/3");
add_line(model,"Eddy/1","HeatOut/1"); add_line(model,"Thermal/1","TempOut/1");
set_param(model,"StopTime","0.3");
simOut=sim(model,"ReturnWorkspaceOutputs","on");
out.heat_density_W_per_m3=simOut.get("native_heat");
out.temperature_K=simOut.get("native_temp");
end

function closeIfLoaded(name)
if bdIsLoaded(name), close_system(name,0); end
end
