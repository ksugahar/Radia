function result = staticElectromagnet(functionName, positional, options)
% Explicit batch entry for radia.static_electromagnet; not a per-step backend.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.static_electromagnet", ...
    functionName, positional, Keywords=options.Keywords);
end
