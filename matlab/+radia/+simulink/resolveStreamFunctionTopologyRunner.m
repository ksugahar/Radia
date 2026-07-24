function runner = resolveStreamFunctionTopologyRunner(modelName,variableName)
%RESOLVESTREAMFUNCTIONTOPOLOGYRUNNER Resolve a checked model/base runner.
arguments
    modelName (1,1) string
    variableName (1,1) string
end
if ~isvarname(variableName)
    error("radia:simulink:StreamFunctionTopologyRunnerVariable", ...
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
        string(runner.Metadata.domain) ~= "stream-function"
    error("radia:simulink:StreamFunctionTopologyRunner", ...
        "Create %s with radia.topopt.makeStreamFunctionAdjointRunner.", ...
        variableName);
end
end
