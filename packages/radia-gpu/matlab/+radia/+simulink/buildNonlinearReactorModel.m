function modelPath = buildNonlinearReactorModel(options)
%BUILDNONLINEARREACTORMODEL Build the public nonlinear toroidal reactor demo.
arguments
    options.ModelName (1,1) string = "radia_nonlinear_reactor"
    options.OutputDirectory (1,1) string = ""
    options.Config (1,1) struct = struct()
    options.Open (1,1) logical = false
end
radia.setup(RequireMex=true);
matlabRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
if strlength(options.OutputDirectory) == 0
    options.OutputDirectory = matlabRoot;
end
if ~isfolder(options.OutputDirectory)
    mkdir(options.OutputDirectory);
end
if bdIsLoaded(options.ModelName)
    close_system(options.ModelName,0);
end
config = options.Config;
if isempty(fieldnames(config))
    config = radia.simulink.makeNonlinearReactorDemoConfig();
end
config = radia.simulink.validateNonlinearReactorConfig(config);

new_system(options.ModelName);
workspace = get_param(options.ModelName,"ModelWorkspace");
workspace.assignin("radia_nonlinear_reactor_config",config);

add_block("simulink/Sources/Sine Wave",options.ModelName+"/Coil_Current", ...
    "Amplitude","20","Frequency","2*pi*50","Bias","0", ...
    "SampleTime",string(config.sample_time_s), ...
    "Position",[35 230 85 280]);
radia.simulink.addNonlinearReactorBlock(options.ModelName, ...
    BlockName="Nonlinear_HDiv_MMM_Reactor",Position=[180 170 500 360]);

add_block("simulink/Signal Routing/Mux",options.ModelName+"/Current_Voltage_Mux", ...
    "Inputs","2","Position",[640 35 645 95]);
add_block("simulink/Sinks/Scope",options.ModelName+"/Current_and_Voltage", ...
    "Position",[760 40 860 90]);
add_block("simulink/Sinks/XY Graph",options.ModelName+"/Flux_Current_Loop", ...
    "Position",[760 125 860 185]);
add_block("simulink/Signal Routing/Mux",options.ModelName+"/Magnetic_State_Mux", ...
    "Inputs","2","Position",[640 220 645 280]);
add_block("simulink/Sinks/Scope",options.ModelName+"/Inductance_and_Peak_B", ...
    "Position",[760 225 860 275]);
radia.simulink.addFieldStatsBlock(options.ModelName, ...
    BlockName="B_Field_Stats",Position=[580 330 720 390]);
add_block("simulink/Sinks/Scope",options.ModelName+"/B_Distribution", ...
    "Position",[760 335 860 385]);
add_block("simulink/Signal Routing/Mux", ...
    options.ModelName+"/Solver_Diagnostics_Mux", ...
    "Inputs","3","Position",[640 430 645 500]);
add_block("simulink/Sinks/Scope",options.ModelName+"/Solver_Diagnostics", ...
    "Position",[760 440 860 490]);

logSignal(connect(options.ModelName,"Coil_Current/1", ...
    "Nonlinear_HDiv_MMM_Reactor/1"),"current_A");
connect(options.ModelName,"Coil_Current/1","Current_Voltage_Mux/1");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/1", ...
    "Current_Voltage_Mux/2"),"voltage_V");
connect(options.ModelName,"Current_Voltage_Mux/1","Current_and_Voltage/1");
connect(options.ModelName,"Coil_Current/1","Flux_Current_Loop/1");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/2", ...
    "Flux_Current_Loop/2"),"flux_Wb_turn");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/3", ...
    "Magnetic_State_Mux/1"),"Ldiff_H");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/4", ...
    "Magnetic_State_Mux/2"),"Bpeak_T");
connect(options.ModelName,"Magnetic_State_Mux/1","Inductance_and_Peak_B/1");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/8", ...
    "B_Field_Stats/1"),"B_samples_T");
connect(options.ModelName,"B_Field_Stats/1","B_Distribution/1");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/5", ...
    "Solver_Diagnostics_Mux/1"),"energy_J");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/6", ...
    "Solver_Diagnostics_Mux/2"),"iterations");
logSignal(connect(options.ModelName,"Nonlinear_HDiv_MMM_Reactor/7", ...
    "Solver_Diagnostics_Mux/3"),"residual");
connect(options.ModelName,"Solver_Diagnostics_Mux/1", ...
    "Solver_Diagnostics/1");

scopeNames = ["Current_and_Voltage","Inductance_and_Peak_B", ...
    "B_Distribution","Solver_Diagnostics"];
for scopeName = scopeNames
    configuration = get_param(options.ModelName+"/"+scopeName, ...
        "ScopeConfiguration");
    configuration.AxesScaling = "Auto";
end
set_param(options.ModelName,"SolverType","Fixed-step", ...
    "Solver","FixedStepDiscrete","FixedStep",string(config.sample_time_s), ...
    "StopTime","0.1","SignalLogging","on", ...
    "SignalLoggingName","logsout","ReturnWorkspaceOutputs","on");
modelPath = fullfile(options.OutputDirectory,options.ModelName+".slx");
save_system(options.ModelName,modelPath);
if options.Open
    open_system(options.ModelName);
else
    close_system(options.ModelName,0);
end
end

function line = connect(modelName,source,destination)
line = add_line(modelName,source,destination,"autorouting","smart");
end

function logSignal(line,name)
sourcePort = get_param(line,"SrcPortHandle");
set_param(line,"Name","");
set_param(sourcePort,"DataLogging","on", ...
    "DataLoggingNameMode","Custom","DataLoggingName",name);
end
