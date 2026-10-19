function modelPath = buildStreamFunctionOptimizationModel(options)
%BUILDSTREAMFUNCTIONOPTIMIZATIONMODEL Build the standalone SF optimizer SLX.
arguments
    options.ModelName (1,1) string = "radia_streamfunction_optimization"
    options.OutputDirectory (1,1) string = ""
    options.ConfigFile (1,1) string = ""
    options.RunRoot (1,1) string = "C:\temp\radia_simulink"
    options.Timeout_s (1,1) double {mustBeNonnegative} = 3600
    options.PythonExecutable (1,1) string = "python"
    options.Runner = []
    options.Solver (1,1) string ...
        {mustBeMember(options.Solver,["mma","sqp"])} = "mma"
    options.SampleTime_s (1,1) double {mustBePositive} = 1
    options.HistoryCapacity (1,1) double ...
        {mustBeInteger,mustBeGreaterThan(options.HistoryCapacity,1)} = 64
    options.Open (1,1) logical = false
end
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
runner = options.Runner;
if isempty(runner)
    runner = radia.topopt.makeStreamFunctionAdjointDemoRunner();
end
if ~isa(runner,"radia.topopt.AdjointRunner") || ...
        ~isfield(runner.Metadata,"domain") || ...
        string(runner.Metadata.domain) ~= "stream-function"
    error("radia:simulink:StreamFunctionTopologyRunner", ...
        "Runner must be created by makeStreamFunctionAdjointRunner.");
end
runner.setSolver(options.Solver);
assignin("base","radia_streamfunction_adjoint_runner",runner);

new_system(options.ModelName);
workspace = get_param(options.ModelName,"ModelWorkspace");
workspace.assignin("radia_streamfunction_adjoint_runner",runner);
blockPath = radia.simulink.addStreamFunctionOptimizationSubsystem( ...
    options.ModelName,"Stream Function Optimization",[190 70 520 190]);
set_param(blockPath, ...
    "config_file",quoteMaskString(options.ConfigFile), ...
    "run_root",quoteMaskString(options.RunRoot), ...
    "timeout_s",compose("%.17g",options.Timeout_s), ...
    "python_executable",quoteMaskString(options.PythonExecutable), ...
    "solver",options.Solver, ...
    "sample_time_s",compose("%.17g",options.SampleTime_s), ...
    "history_capacity",string(options.HistoryCapacity));

add_block("simulink/Sources/Constant",options.ModelName + "/Run Design", ...
    Value="false",OutDataTypeStr="boolean",Position=[35 85 95 115]);
add_block("simulink/Sources/Constant",options.ModelName + "/Run Adjoint", ...
    Value="false",OutDataTypeStr="boolean",Position=[35 145 95 175]);
add_line(options.ModelName,"Run Design/1", ...
    "Stream Function Optimization/1");
add_line(options.ModelName,"Run Adjoint/1", ...
    "Stream Function Optimization/2");

outputNames = [ ...
    "design_status","primary","elapsed_s", ...
    "objective","optimization_status","iterations","design", ...
    "history_count","history_iteration","history_objective", ...
    "history_constraint"];
for index = 1:numel(outputNames)
    y = 25 + (index-1)*35;
    add_block("simulink/Ports & Subsystems/Out1", ...
        options.ModelName + "/" + outputNames(index),Port=string(index), ...
        Position=[690 y 720 y+20]);
    add_line(options.ModelName, ...
        "Stream Function Optimization/" + index, ...
        outputNames(index) + "/1");
end
set_param(options.ModelName,"Solver","FixedStepDiscrete", ...
    "FixedStep",compose("%.17g",options.SampleTime_s),"StopTime","1");
modelPath = fullfile(options.OutputDirectory,options.ModelName + ".slx");
save_system(options.ModelName,modelPath);
if options.Open
    open_system(options.ModelName);
else
    close_system(options.ModelName,0);
end
end

function value = quoteMaskString(value)
value = replace(string(value),"'","''");
value = "'" + value + "'";
end
