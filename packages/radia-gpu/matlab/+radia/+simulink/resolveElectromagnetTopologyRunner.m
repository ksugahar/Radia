function runner = resolveElectromagnetTopologyRunner(modelName,variableName)
%RESOLVEELECTROMAGNETTOPOLOGYRUNNER Resolve a checked model/base runner.
arguments
    modelName (1,1) string
    variableName (1,1) string
end
if ~isvarname(variableName)
    error("radia:simulink:ElectromagnetTopologyRunnerVariable", ...
        "Runner must be named by a valid MATLAB variable.");
end
runner = [];
if bdIsLoaded(modelName)
    workspace = get_param(modelName,"ModelWorkspace");
    if workspace.hasVariable(variableName)
        runner = workspace.getVariable(variableName);
    end
end
if isempty(runner) && ...
        evalin("base","exist('" + variableName + "','var')") == 1
    runner = evalin("base",variableName);
end
if ~isa(runner,"radia.topopt.AdjointRunner") || ...
        ~isfield(runner.Metadata,"domain") || ...
        string(runner.Metadata.domain) ~= "electromagnet-topology"
    error("radia:simulink:ElectromagnetTopologyRunner", ...
        "Create %s with radia.topopt.makeElectromagnetAdjointRunner.", ...
        variableName);
end
end
