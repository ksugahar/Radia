function result = analytical(moduleName, functionName, positional, options)
%ANALYTICAL Call a Python-only analytic formula module.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
if ~startsWith(moduleName, "analytical_formulas.") && ...
        ~ismember(moduleName, ["analytical_magnet","cylindrical_magnet","round_bodies"])
    error("radia:python:Module", "Unsupported analytical module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
