function result = electromagnetics(moduleName, functionName, positional, options)
%ELECTROMAGNETICS Call another Python-only electromagnetic method operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = any(startsWith(moduleName, ["planar_","vim."])) || ...
    ismember(moduleName, ["axifem","biot_savart","clebsch_potential", ...
    "hysteresis_io","radia_ngsolve","soft_iron","veriloga_generator", ...
    "vim"]);
if ~allowed
    error("radia:python:Module", ...
        "Unsupported electromagnetic module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
