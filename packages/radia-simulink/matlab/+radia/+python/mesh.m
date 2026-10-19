function result = mesh(moduleName, functionName, positional, options)
%MESH Call a Python-only CAD, mesh import, or export operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = contains(moduleName, ["mesh","cubit","gmsh","step"]);
if ~any(allowed)
    error("radia:python:Module", "Unsupported CAD/mesh module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
