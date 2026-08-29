function modelPath = buildOptunaTeachingModel(options)
%BUILDOPTUNATEACHINGMODEL Build the standalone student Optuna laboratory.
%   The same model layout supports a known optimum, a biobjective Pareto
%   exercise, and a deterministic complete/pruned/failed exercise.
arguments
    options.OutputPath (1,1) string = ""
    options.Exercise (1,1) string {mustBeMember(options.Exercise, ...
        ["quadratic","pareto","reliability"])} = "quadratic"
    options.StoragePath (1,1) string = ""
    options.Save (1,1) logical = true
end

matlabRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
if strlength(options.OutputPath) == 0
    options.OutputPath = fullfile(matlabRoot, ...
        "radia_optuna_teaching.slx");
end
[folder, modelName, extension] = fileparts(options.OutputPath);
if extension == ""
    options.OutputPath = options.OutputPath + ".slx";
elseif ~strcmpi(extension, ".slx")
    error("radia:simulink:TeachingModelPath", ...
        "OutputPath must name an SLX file.");
end
if strlength(folder) > 0 && ~isfolder(folder)
    mkdir(folder);
end
modelName = string(matlab.lang.makeValidName(modelName));
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
new_system(modelName);
cleanup = onCleanup(@()closeIfLoaded(modelName));

[objectiveName, directions, parameters, sampler, pruner, nTrials] = ...
    exerciseConfiguration(options.Exercise);
workspace = get_param(modelName, "ModelWorkspace");
assignin(workspace, "radia_optuna_teaching_parameters", parameters);
assignin(workspace, "radia_optuna_teaching_exercise", options.Exercise);
storagePath = options.StoragePath;
if strlength(storagePath) == 0
    storagePath = "C:\temp\" + modelName + "_" + ...
        options.Exercise + ".mat";
end

blockPath = radia.simulink.buildOptunaBlock(modelName, ...
    ObjectiveFcn=objectiveName, NumTrials=nTrials, ...
    Directions=directions, StoragePath=storagePath, ...
    SampleTime_s=0.1, Sampler=sampler, Seed=20260829, ...
    Pruner=pruner, ...
    ParameterSpec="radia_optuna_teaching_parameters", ...
    ModelName=modelName, LiveVisualization=false, Save=false);
set_param(blockPath, "Position", [300 90 600 430]);

inputNames = ["Start","Cancel","Pause","Resume", ...
    "Selected Trial","Apply"];
inputValues = ["1","0","0","0","0","0"];
for index = 1:numel(inputNames)
    path = modelName + "/" + inputNames(index);
    add_block("simulink/Sources/Constant", path, ...
        "Value", inputValues(index), ...
        "Position", [55, 80 + 58*index, 145, 110 + 58*index]);
    add_line(modelName, inputNames(index) + "/1", ...
        "Optuna Optimization/" + index);
end

workspaceOutputs = struct( ...
    "best", 1, "status", 3, "completed", 4, ...
    "pareto_count", 8, "failed", 12, "attempted", 13, ...
    "selected", 15, "pruned", 16, "checkpoint", 18);
fields = fieldnames(workspaceOutputs);
connected = false(1, 18);
for index = 1:numel(fields)
    field = string(fields{index});
    port = workspaceOutputs.(fields{index});
    sink = modelName + "/Log " + field;
    add_block("simulink/Sinks/To Workspace", sink, ...
        "VariableName", "teaching_" + field, "SaveFormat", "Array", ...
        "Position", [720, 45 + 42*index, 845, 70 + 42*index]);
    add_line(modelName, "Optuna Optimization/" + port, ...
        "Log " + field + "/1");
    connected(port) = true;
end
for port = find(~connected)
    sink = modelName + "/Unused " + port;
    add_block("simulink/Sinks/Terminator", sink, ...
        "Position", [900, 35 + 24*port, 920, 55 + 24*port]);
    add_line(modelName, "Optuna Optimization/" + port, ...
        "Unused " + port + "/1");
end

displayPorts = [1, 3, 13, 16];
displayNames = ["Best Value","Session Status","Attempted Trials", ...
    "Pruned Trials"];
for index = 1:numel(displayPorts)
    path = modelName + "/" + displayNames(index);
    add_block("simulink/Sinks/Display", path, ...
        "Position", [980, 80 + 75*index, 1090, 115 + 75*index]);
    add_line(modelName, "Optuna Optimization/" + displayPorts(index), ...
        displayNames(index) + "/1");
end

set_param(modelName, "Solver", "FixedStepDiscrete", ...
    "FixedStep", "0.1", "StopTime", compose("%.1f", 0.1*(nTrials+2)), ...
    "ModelBrowserVisibility", "on");
annotation = Simulink.Annotation(modelName, ...
    "Radia Optuna teaching lab: " + options.Exercise + newline + ...
    "Run, inspect the trial table, change sampler/seed/search space, " + ...
    "then select and apply a completed trial.");
annotation.Position = [300, 25, 790, 65];
annotation.FontSize = 11;

if options.Save
    save_system(modelName, options.OutputPath);
end
modelPath = options.OutputPath;
clear cleanup
end

function [objectiveName, directions, parameters, sampler, pruner, nTrials] = ...
        exerciseConfiguration(exercise)
switch exercise
    case "quadratic"
        objectiveName = "radia.optuna.teachingQuadraticObjective";
        directions = "minimize";
        parameters = radia.optuna.OptimizationParameter("x", ...
            Value=0, Minimum=-1, Maximum=1);
        sampler = "tpe";
        pruner = "median";
        nTrials = 12;
    case "pareto"
        objectiveName = "radia.optuna.teachingParetoObjective";
        directions = ["minimize","minimize"];
        parameters = radia.optuna.OptimizationParameter("x", ...
            Value=0.5, Minimum=0, Maximum=1);
        sampler = "nsgaii";
        pruner = "none";
        nTrials = 12;
    otherwise
        objectiveName = "radia.optuna.teachingReliabilityObjective";
        directions = "minimize";
        parameters = radia.optuna.OptimizationParameter("x", ...
            Value=0, Minimum=-1, Maximum=1, Step=0.5);
        sampler = "bruteforce";
        pruner = "none";
        nTrials = 5;
end
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
