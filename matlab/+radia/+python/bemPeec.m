function result = bemPeec(moduleName, functionName, positional, options)
%BEMPEEC Call a high-level Python BEM, PEEC, SIBC, or dielectric operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = any(startsWith(moduleName, ["bem","peec","ngsbem_"])) || ...
    ismember(moduleName, ["dielectric_solver","fasthenry_parser", ...
    "scalar_bie_sibc"]);
if ~allowed
    error("radia:python:Module", "Unsupported BEM/PEEC module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
