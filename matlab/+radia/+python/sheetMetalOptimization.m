function result = sheetMetalOptimization(functionName, positional, options)
%SHEETMETALOPTIMIZATION Call the Python NGSolve/Cubit shape API explicitly.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.sheet_metal_optimization", ...
    functionName, positional, Keywords=options.Keywords);
end
