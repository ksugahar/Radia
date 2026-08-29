function parameters = getParameterFromModel(model, names)
%GETPARAMETERFROMMODEL Read tunable parameters out of a Simulink model.
%   PARAMETERS = RADIA.OPTUNA.GETPARAMETERFROMMODEL(MODEL, NAMES) returns a
%   radia.optuna.OptimizationParameter for each name, with Value taken from
%   the model, matching how sdo.getParameterFromModel starts an optimization.
%
%   Bounds are deliberately left at -Inf/Inf. sdo does the same, and
%   radia.optuna.optimize refuses to search a parameter whose bounds are not
%   finite: a study needs a range, and inventing one from the current value
%   would quietly decide the search space on the user's behalf.
%
%   Each name is resolved in the model workspace first, then the base
%   workspace, matching how Simulink resolves a block parameter expression. A
%   name that resolves in neither is an error naming both places looked, not
%   a parameter silently dropped from the returned array.
%
%   Example:
%       p = radia.optuna.getParameterFromModel("pid_plant", ["Kp" "Ki"]);
%       p(1).Minimum = 0;   p(1).Maximum = 10;
%       p(2).Minimum = 1e-3; p(2).Maximum = 1; p(2).Transform = "log";
%       [x, fval] = radia.optuna.optimize(@cost, p, opts);
%
%   See also radia.optuna.optimize, radia.optuna.OptimizationParameter.

arguments
    model (1,1) string
    names (1,:) string
end
if exist("get_param", "file") ~= 2
    error("radia:optuna:SimulinkRequired", ...
        "getParameterFromModel needs Simulink; get_param is not available.");
end

wasLoaded = bdIsLoaded(model);
if ~wasLoaded
    load_system(model);
end
cleanup = onCleanup(@() closeIfOpened(model, wasLoaded));

workspace = get_param(model, "ModelWorkspace");
parameters = radia.optuna.OptimizationParameter.empty(1, 0);
for index = 1:numel(names)
    name = names(index);
    value = [];
    found = false;
    if ~isempty(workspace) && workspace.hasVariable(char(name))
        value = workspace.getVariable(char(name));
        found = true;
    elseif evalin("base", sprintf("exist('%s','var')", name)) == 1
        value = evalin("base", char(name));
        found = true;
    end
    if ~found
        error("radia:optuna:ParameterNotFound", ...
            "Parameter '%s' is in neither the '%s' model workspace nor " + ...
            "the base workspace.", name, model);
    end
    if ~isnumeric(value) || ~isscalar(value)
        error("radia:optuna:ParameterNotScalar", ...
            "Parameter '%s' is %s of size %s; optimize searches scalar " + ...
            "numeric parameters.", name, class(value), ...
            mat2str(size(value)));
    end
    parameters(index) = radia.optuna.OptimizationParameter(name, ...
        Value=double(value));
end
end

function closeIfOpened(model, wasLoaded)
if ~wasLoaded && bdIsLoaded(model)
    close_system(model, 0);
end
end
