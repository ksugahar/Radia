function result = electromagnetValidation(functionName, positional, options)
% Explicit batch entry for radia.electromagnet_validation; not a per-step backend.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.electromagnet_validation", ...
    functionName, positional, Keywords=options.Keywords);
end
