function result = inductionHeating(moduleName, functionName, positional, options)
%INDUCTIONHEATING Call a Python-only ESIM or IH batch operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
if ~any(startsWith(moduleName, ["esim_","ih_"])) && ...
        ~ismember(moduleName, ["em_design","em_material"])
    error("radia:python:Module", "Unsupported ESIM/IH module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
