function files = buildLTspicePIDOptunaExample(options)
%BUILDLTSPICEPIDOPTUNAEXAMPLE Build the LTspice PID Optuna Simulink example.
arguments
    options.OutputDirectory (1,1) string = radia.simulink.exampleDirectory()
    options.NumTrials (1,1) double {mustBeInteger,mustBePositive} = 8
    options.OpenModel (1,1) logical = false
end

outputDirectory = string(options.OutputDirectory);
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
plantName = "radia_ltspice_pid_plant";
harnessName = "radia_ltspice_pid_optuna";
closeIfLoaded([plantName,harnessName]);
buildPlant(plantName,outputDirectory);
plantFile=fullfile(outputDirectory,plantName+".slx");
buildHarness(harnessName,outputDirectory,options.NumTrials);
files = struct("harness",fullfile(outputDirectory,harnessName+".slx"), ...
    "plant",plantFile);
if options.OpenModel, open_system(harnessName); end
end

function buildPlant(model,outputDirectory)
new_system(model); cleanup = onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Step",model+"/Reference", ...
    Time="0",Before="0",After="1",SampleTime="0.005");
add_block("simulink/Math Operations/Sum",model+"/Error",Inputs="+-");
add_block("simulink/Discrete/Discrete PID Controller",model+"/PID", ...
    Controller="PID",P="Kp",I="Ki",D="Kd",N="1000", ...
    SampleTime="0.005");
add_block("simulink/Discontinuities/Saturation",model+"/Command limit", ...
    UpperLimit="5",LowerLimit="0");
netlist = fullfile(radia.simulink.exampleDirectory(), ...
    "samples","ltspice_pid_rc_plant.cir");
if ~isfile(netlist)
    error("radia:simulink:LTspicePIDNetlist", ...
        "Bundled LTspice PID plant netlist was not found: %s",netlist);
end
radia.simulink.buildLTspiceBlock(model,Netlist=netlist, ...
    InputNames="control",OutputTraces="V(output)", ...
    SampleTime_s=5e-3,MaxStep_s=1e-4,Timeout_s=30,Save=false);
add_block("simulink/Sinks/To Workspace",model+"/Plant output", ...
    VariableName="pid_output",SaveFormat="Timeseries");
add_block("simulink/Sinks/To Workspace",model+"/Tracking error", ...
    VariableName="pid_error",SaveFormat="Timeseries");
add_block("simulink/Sinks/To Workspace",model+"/Control effort", ...
    VariableName="pid_control",SaveFormat="Timeseries");
add_block("simulink/Discrete/Unit Delay",model+"/Feedback state", ...
    InitialCondition="0",SampleTime="0.005");
add_block("simulink/Signal Routing/Mux",model+"/Scope Mux",Inputs="3");
add_block("simulink/Sinks/Scope",model+"/Closed-loop scope");
add_line(model,"Reference/1","Error/1");
add_line(model,"Error/1","PID/1"); add_line(model,"Error/1","Tracking error/1");
add_line(model,"PID/1","Command limit/1");
add_line(model,"Command limit/1","LTspice Circuit/1");
add_line(model,"Command limit/1","Control effort/1");
add_line(model,"LTspice Circuit/1","Feedback state/1");
add_line(model,"Feedback state/1","Error/2");
add_line(model,"LTspice Circuit/1","Plant output/1");
add_line(model,"Reference/1","Scope Mux/1");
add_line(model,"LTspice Circuit/1","Scope Mux/2");
add_line(model,"Command limit/1","Scope Mux/3");
add_line(model,"Scope Mux/1","Closed-loop scope/1");
set_param(model,Solver="FixedStepDiscrete",FixedStep="0.005",StopTime="0.025");
set_param(model,"PreLoadFcn","Kp=1; Ki=50; Kd=1e-4;");
Simulink.BlockDiagram.arrangeSystem(model);
save_system(model,fullfile(outputDirectory,model+".slx"));
clear cleanup; close_system(model,0);
end

function buildHarness(model,outputDirectory,numTrials)
new_system(model); cleanup = onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Step",model+"/Start optimization", ...
    Time="0",Before="0",After="1",SampleTime="0.1");
add_block("simulink/Sources/Constant",model+"/Cancel optimization", ...
    Value="0",SampleTime="0.1");
radia.simulink.buildOptunaStudyBlock(model, ...
    BlockName="Optuna Study", ...
    ObjectiveFcn="radia.simulink.ltspicePIDObjective", ...
    NumTrials=numTrials,SampleTime_s=0.1,Sampler="cmaes",Save=false);
add_block("simulink/Sinks/Display",model+"/Best objective");
add_block("simulink/Sinks/Display",model+"/Completed trials");
add_block("simulink/Sinks/Display",model+"/Status");
add_block("simulink/Sinks/Display",model+"/Best trial");
add_line(model,"Start optimization/1","Optuna Study/1");
add_line(model,"Cancel optimization/1","Optuna Study/2");
add_line(model,"Optuna Study/1","Best objective/1");
add_line(model,"Optuna Study/2","Status/1");
add_line(model,"Optuna Study/3","Completed trials/1");
add_line(model,"Optuna Study/4","Best trial/1");
set_param(model,Solver="FixedStepDiscrete",FixedStep="0.1", ...
    StopTime=compose("%.17g",numTrials*0.1));
Simulink.BlockDiagram.arrangeSystem(model);
save_system(model,fullfile(outputDirectory,model+".slx"));
clear cleanup; close_system(model,0);
end

function closeIfLoaded(models)
for model = reshape(string(models),1,[])
    if bdIsLoaded(model), close_system(model,0); end
end
end
