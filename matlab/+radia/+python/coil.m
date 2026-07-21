function result = coil(moduleName, functionName, positional, options)
%COIL Call a Python-only coil construction or topology operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
if ~startsWith(moduleName, "coil") && ...
        ~ismember(moduleName, ["coils","filament_bundle","workpiece_surface"])
    error("radia:python:Module", "Unsupported coil module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
