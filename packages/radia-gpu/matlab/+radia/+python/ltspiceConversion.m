function result = ltspiceConversion(moduleName, functionName, positional, options)
%LTSPICECONVERSION Call a Radia LTspice conversion or topology operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end

allowedModules = ["", "conversion", "topology", "parser", ...
    "parser.asc_parser", "parser.asc_to_schemdraw", ...
    "parser.cir_to_schemdraw", "parser.netlist_to_asc", ...
    "parser.schemdraw_to_cir", "parser.schemdraw_to_ltspice"];
if ~ismember(moduleName, allowedModules)
    error("radia:python:Module", ...
        "Unsupported Radia LTspice conversion module: %s", moduleName);
end

pythonModule = "radia.ltspice";
if strlength(moduleName) > 0
    pythonModule = pythonModule + "." + moduleName;
end
result = radia.internal.callPython(pythonModule, functionName, positional, ...
    Keywords=options.Keywords);
end
