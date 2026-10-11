function [runner, study] = resolveStreamFunctionOptunaObjects( ...
        modelName, runnerVariable, studyVariable)
%RESOLVESTREAMFUNCTIONOPTUNAOBJECTS Resolve checked model/base objects.
arguments
    modelName (1,1) string
    runnerVariable (1,1) string
    studyVariable (1,1) string
end
if ~isvarname(runnerVariable) || ~isvarname(studyVariable)
    error("radia:simulink:StreamFunctionOptunaVariable", ...
        "Runner and study must be named by valid MATLAB variables.");
end
runner = resolveVariable(modelName, runnerVariable);
study = resolveVariable(modelName, studyVariable);
if ~isa(runner, "radia.stream.OptunaRunner")
    error("radia:simulink:StreamFunctionOptunaRunner", ...
        "Create %s with radia.stream.OptunaRunner.", runnerVariable);
end
if ~isa(study, "radia.optuna.Study")
    error("radia:simulink:StreamFunctionOptunaStudy", ...
        "Create %s with radia.optuna.createStudy.", studyVariable);
end
end

function value = resolveVariable(modelName, variableName)
value = [];
if bdIsLoaded(modelName)
    workspace = get_param(modelName, "ModelWorkspace");
    if workspace.hasVariable(variableName)
        value = workspace.getVariable(variableName);
    end
end
if isempty(value) && ...
        evalin("base", "exist('" + variableName + "','var')") == 1
    value = evalin("base", variableName);
end
end
