function result = optimization(moduleName, functionName, positional, options)
%OPTIMIZATION Call a Python-only topology or application optimization operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = contains(moduleName, "optimization") || ...
    ismember(moduleName, ["hcurl_topology_optimization","ih_optimize"]);
if ~allowed
    error("radia:python:Module", "Unsupported optimization module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
