function result = esrfExamples(functionName, positional, options)
% Explicit batch entry for radia.esrf_examples; not a per-step backend.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.esrf_examples", ...
    functionName, positional, Keywords=options.Keywords);
end
