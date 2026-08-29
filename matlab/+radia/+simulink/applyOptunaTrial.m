function values = applyOptunaTrial(blockPath, trialNumber)
%APPLYOPTUNATRIAL Apply a saved completed trial without Simulink wiring.
arguments
    blockPath (1,1) string
    trialNumber (1,1) double {mustBeInteger}
end
storagePath = evaluatedMaskValue(blockPath, "storage_path");
if strlength(storagePath) == 0
    error("radia:simulink:OptunaStudyStorage", ...
        "Set Study MAT file before applying a trial.");
end
study = radia.optuna.loadStudy(storage=storagePath);
if trialNumber == -1
    best = study.bestSolution();
    if ~best.available
        error("radia:simulink:OptunaStudySelection", ...
            "The study has no completed feasible trial to apply.");
    end
    trialNumber = best.trial_number;
elseif trialNumber < 0
    error("radia:simulink:OptunaStudySelection", ...
        "Trial must be nonnegative, or -1 for the best trial.");
end
sessionPath = storagePath + ".session.mat";
modelName = evaluatedMaskValue(blockPath, "model_name");
if strlength(modelName) == 0
    modelName = string(bdroot(blockPath));
end
if isfile(sessionPath)
    session = radia.optuna.OptimizationSession.load( ...
        sessionPath, ModelName=modelName);
    session.selectTrial(trialNumber);
    values = session.applySelectedToModel( ...
        ModelName=modelName, Target="model");
else
    completed = study.trials("COMPLETE");
    row = find(completed.TrialNumber == trialNumber, 1);
    if isempty(row)
        error("radia:simulink:OptunaStudySelection", ...
            "Trial %d is not a completed trial.", trialNumber);
    end
    values = completed.Params{row};
    applyParameterValues(values, modelName);
end
set_param(blockPath, "selected_trial", compose("%.0f", trialNumber));
fprintf("Applied Optuna trial %d to model '%s'.\n", ...
    trialNumber, modelName);
end

function applyParameterValues(values, modelName)
if ~isstruct(values) || ~isscalar(values)
    error("radia:simulink:OptunaStudySelection", ...
        "The selected trial has no applicable parameter struct.");
end
if ~bdIsLoaded(modelName)
    load_system(modelName);
end
workspace = get_param(modelName, "ModelWorkspace");
names = string(fieldnames(values));
for index = 1:numel(names)
    name = names(index);
    value = values.(name);
    if workspace.hasVariable(char(name))
        workspace.assignin(char(name), value);
    else
        assignin("base", char(name), value);
    end
end
end

function value = evaluatedMaskValue(blockPath, name)
variables = get_param(blockPath, "MaskWSVariables");
match = find(strcmp({variables.Name}, name), 1);
if isempty(match)
    error("radia:simulink:OptunaStudyMask", ...
        "Block '%s' has no '%s' mask parameter.", blockPath, name);
end
value = string(variables(match).Value);
end
