function result = openBoundary(moduleName, functionName, positional, options)
%OPENBOUNDARY Call a Python-only Kelvin, DtN, or potential operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = any(startsWith(moduleName, ["kelvin_","open_boundary."])) || ...
    ismember(moduleName, ["cohomology","cohomology_cut","ima_field", ...
    "infinite_element","scalar_potential_solver","vector_potential_solver"]);
if ~allowed
    error("radia:python:Module", "Unsupported open-boundary module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
