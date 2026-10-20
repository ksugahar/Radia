function modelFile = buildPEECCoilRLCModel(options)
%BUILDPEECCOILRLCMODEL Build an explicit PEEC R-L coil with fixed C.
arguments
    options.OutputDirectory (1,1) string = radia.simulink.exampleDirectory()
    options.ModelName (1,1) string = "radia_peec_coil_rlc_plant"
    options.OpenModel (1,1) logical = false
end

if ~isfolder(options.OutputDirectory)
    mkdir(options.OutputDirectory);
end
model = options.ModelName;
modelFile = fullfile(options.OutputDirectory,model+".slx");
if bdIsLoaded(model)
    close_system(model,0);
end
if isfile(modelFile)
    delete(modelFile);
end

new_system(model);
cleanup = onCleanup(@() closeIfLoaded(model));
set_param(model, ...
    "Solver","ode45", ...
    "StopTime","0.012", ...
    "MaxStep","1e-5", ...
    "RelTol","1e-9", ...
    "AbsTol","1e-11", ...
    "ReturnWorkspaceOutputs","on", ...
    "PreLoadFcn","R_ohm=0.02; L_H=8e-6; C_F=3.3e-3;", ...
    "ModelBrowserVisibility","off");

add_block("simulink/Sources/Constant",model+"/Zero source voltage", ...
    "Value","0","Position",[35 145 95 175]);
add_block("simulink/Math Operations/Sum",model+"/Voltage balance", ...
    "Inputs","+--","Position",[150 120 180 200]);
add_block("simulink/Math Operations/Gain", ...
    model+"/PEEC Inductance 1 over L", ...
    "Gain","1/L_H","Position",[235 135 365 185]);
add_block("simulink/Continuous/Integrator",model+"/Coil current", ...
    "InitialCondition","0","Position",[415 140 455 180]);
add_block("simulink/Math Operations/Gain",model+"/Copper resistance R", ...
    "Gain","R_ohm","Position",[515 235 630 275]);
add_block("simulink/Math Operations/Gain",model+"/Fixed capacitance 1 over C", ...
    "Gain","1/C_F","Position",[515 75 650 115]);
add_block("simulink/Continuous/Integrator",model+"/Capacitor voltage", ...
    "InitialCondition","1","Position",[700 75 740 115]);

add_block("simulink/Sinks/To Workspace",model+"/Ring voltage log", ...
    "VariableName","ring_voltage","SaveFormat","Timeseries", ...
    "Position",[805 55 925 85]);
add_block("simulink/Sinks/To Workspace",model+"/Coil current log", ...
    "VariableName","coil_current","SaveFormat","Timeseries", ...
    "Position",[515 145 635 175]);
add_block("simulink/Ports & Subsystems/Out1", ...
    model+"/capacitor_voltage_V","Port","1", ...
    "Position",[805 95 925 115]);
add_block("simulink/Ports & Subsystems/Out1", ...
    model+"/coil_current_A","Port","2", ...
    "Position",[515 185 635 205]);

add_block("simulink/Sources/Constant",model+"/PEEC inductance H", ...
    "Value","L_H","Position",[45 330 150 360]);
add_block("simulink/Sinks/To Workspace",model+"/PEEC inductance log", ...
    "VariableName","peec_inductance","SaveFormat","Timeseries", ...
    "Position",[215 325 350 365]);
add_block("simulink/Sources/Constant",model+"/PEEC resistance Ohm", ...
    "Value","R_ohm","Position",[45 390 150 420]);
add_block("simulink/Sinks/To Workspace",model+"/PEEC resistance log", ...
    "VariableName","peec_resistance","SaveFormat","Timeseries", ...
    "Position",[215 385 350 425]);

add_block("simulink/Signal Routing/Mux",model+"/Scope signals", ...
    "Inputs","2","Position",[780 180 785 240]);
add_block("simulink/Sinks/Scope",model+"/RLC scope", ...
    "Position",[850 190 885 225]);

add_line(model,"Zero source voltage/1","Voltage balance/1");
add_line(model,"Voltage balance/1","PEEC Inductance 1 over L/1");
add_line(model,"PEEC Inductance 1 over L/1","Coil current/1");
add_line(model,"Coil current/1","Copper resistance R/1");
add_line(model,"Copper resistance R/1","Voltage balance/2");
add_line(model,"Coil current/1","Fixed capacitance 1 over C/1");
add_line(model,"Fixed capacitance 1 over C/1","Capacitor voltage/1");
add_line(model,"Capacitor voltage/1","Voltage balance/3");
add_line(model,"Capacitor voltage/1","Ring voltage log/1");
add_line(model,"Capacitor voltage/1","capacitor_voltage_V/1");
add_line(model,"Coil current/1","Coil current log/1");
add_line(model,"Coil current/1","coil_current_A/1");
add_line(model,"PEEC inductance H/1","PEEC inductance log/1");
add_line(model,"PEEC resistance Ohm/1","PEEC resistance log/1");
add_line(model,"Capacitor voltage/1","Scope signals/1");
add_line(model,"Coil current/1","Scope signals/2");
add_line(model,"Scope signals/1","RLC scope/1");

annotationText = "Geometry is evaluated outside the time-step loop." + newline + ...
    "Radia HACApK PEEC supplies L_H and copper R_ohm; C_F stays fixed.";
annotation = Simulink.Annotation(model,char(annotationText));
annotation.Position = [35 25 650 65];
Simulink.BlockDiagram.arrangeSystem(model);
save_system(model,modelFile);
if options.OpenModel
    clear cleanup
    open_system(model);
else
    close_system(model,0);
    clear cleanup
end
end

function closeIfLoaded(model)
if bdIsLoaded(model)
    close_system(model,0);
end
end
