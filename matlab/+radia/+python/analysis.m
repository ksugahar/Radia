function result = analysis(moduleName, functionName, positional, options)
%ANALYSIS Call a Python-only analysis or model-reduction operation.
arguments
    moduleName (1,1) string {mustBeMember(moduleName, ...
        ["analysis","lanczos_reduction","prima_hacapk"])}
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
