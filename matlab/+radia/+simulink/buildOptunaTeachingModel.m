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

[objectiveName, directions, parameterExpression, sampler, pruner, ...
    nTrials, initialValue] = ...
    exerciseConfiguration(options.Exercise);
workspace = get_param(modelName, "ModelWorkspace");
assignin(workspace, "radia_optuna_teaching_exercise", options.Exercise);
assignin(workspace, "x", initialValue);
storagePath = options.StoragePath;
if strlength(storagePath) == 0
    storagePath = "C:\temp\" + modelName + "_" + ...
        options.Exercise + ".mat";
end

blockPath = radia.simulink.buildOptunaStudyBlock(modelName, ...
    BlockName="Optuna Study", ...
    ObjectiveFcn=objectiveName, NumTrials=nTrials, ...
    Directions=directions, StoragePath=storagePath, ...
    SampleTime_s=0.1, Sampler=sampler, Seed=20260829, ...
    Pruner=pruner, ...
    ParameterSpec=parameterExpression, ...
    ModelName=modelName, LiveVisualization=false, Save=false);
set_param(blockPath, "Position", [300 120 560 300]);

inputNames = ["Start","Cancel"];
inputValues = ["1","0"];
for index = 1:numel(inputNames)
    path = modelName + "/" + inputNames(index);
    add_block("simulink/Sources/Constant", path, ...
        "Value", inputValues(index), ...
        "Position", [80, 105 + 90*index, 170, 135 + 90*index]);
    add_line(modelName, inputNames(index) + "/1", ...
        "Optuna Study/" + index);
end

dashboard = addTeachingDashboard(modelName);
for port = 1:4
    add_line(modelName, "Optuna Study/" + port, ...
        "Study Dashboard/" + port);
end

set_param(modelName, "Solver", "FixedStepDiscrete", ...
    "FixedStep", "0.1", "StopTime", compose("%.1f", 0.1*(nTrials+2)), ...
    "ModelBrowserVisibility", "on");
annotation = Simulink.Annotation(modelName, ...
    "Radia Optuna teaching lab: " + options.Exercise + newline + ...
    "Run, review the saved study from the block mask, change sampler, " + ...
    "seed, pruner, bounds, or budget, then run again without rewiring.");
annotation.Position = [275, 35, 925, 75];
annotation.FontSize = 11;

if options.Save
    save_system(modelName, options.OutputPath);
end
modelPath = options.OutputPath;
clear cleanup
end

function dashboard = addTeachingDashboard(modelName)
dashboard = modelName + "/Study Dashboard";
add_block("simulink/Ports & Subsystems/Subsystem", dashboard, ...
    Position=[700 120 920 300]);
delete_line(dashboard, "In1/1", "Out1/1");
delete_block(dashboard + "/In1");
delete_block(dashboard + "/Out1");
inputNames = ["Best","Status","Progress","Best Trial"];
variableNames = ["best","status","attempted","best_trial"];
displayNames = ["Best Value","Session Status","Attempted Trials", ...
    "Best Trial Number"];
for index = 1:4
    add_block("simulink/Ports & Subsystems/In1", ...
        dashboard + "/" + inputNames(index), Port=string(index), ...
        Position=[30, 45 + 70*index, 60, 59 + 70*index]);
    add_block("simulink/Sinks/To Workspace", ...
        dashboard + "/Log " + variableNames(index), ...
        VariableName="teaching_" + variableNames(index), ...
        SaveFormat="Array", ...
        Position=[145, 30 + 70*index, 270, 55 + 70*index]);
    add_block("simulink/Sinks/Display", ...
        dashboard + "/" + displayNames(index), ...
        Position=[325, 55 + 70*index, 455, 90 + 70*index]);
    add_line(dashboard, inputNames(index) + "/1", ...
        "Log " + variableNames(index) + "/1");
    add_line(dashboard, inputNames(index) + "/1", ...
        displayNames(index) + "/1");
end
mask = Simulink.Mask.create(dashboard);
mask.Description = "Compact live progress plus workspace logs. Full trial " + ...
    "history remains in the Optuna Study MAT file.";
mask.Display = "disp('Study Dashboard');" + ...
    "port_label('input',1,'best');port_label('input',2,'status');" + ...
    "port_label('input',3,'progress');" + ...
    "port_label('input',4,'best trial');";
set_param(dashboard, UserData=struct( ...
    "role", "student-progress-dashboard", ...
    "full_history", "normalized-mat-study", ...
    "top_level_output_lines", 4), UserDataPersistent="on");
end

function [objectiveName, directions, parameterExpression, sampler, pruner, ...
        nTrials, initialValue] = ...
        exerciseConfiguration(exercise)
switch exercise
    case "quadratic"
        objectiveName = "radia.optuna.teachingQuadraticObjective";
        directions = "minimize";
        parameterExpression = ...
            "radia.optuna.OptimizationParameter('x'," + ...
            "Value=0,Minimum=-1,Maximum=1)";
        sampler = "tpe";
        pruner = "median";
        nTrials = 12;
        initialValue = 0;
    case "pareto"
        objectiveName = "radia.optuna.teachingParetoObjective";
        directions = ["minimize","minimize"];
        parameterExpression = ...
            "radia.optuna.OptimizationParameter('x'," + ...
            "Value=0.5,Minimum=0,Maximum=1)";
        sampler = "nsgaii";
        pruner = "none";
        nTrials = 12;
        initialValue = 0.5;
    otherwise
        objectiveName = "radia.optuna.teachingReliabilityObjective";
        directions = "minimize";
        parameterExpression = ...
            "radia.optuna.OptimizationParameter('x'," + ...
            "Value=0,Minimum=-1,Maximum=1,Step=0.5)";
        sampler = "bruteforce";
        pruner = "none";
        nTrials = 5;
        initialValue = 0;
end
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
