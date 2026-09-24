function result = meshedCurrent(functionName, positional, options)
% Explicit batch entry for radia.meshed_current; not a per-step backend.
% Python NGSolve objects remain in result.value; native MEX handles cannot
% be passed to this interface. Parallel assembly uses a caller-owned TaskManager.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.meshed_current", ...
    functionName, positional, Keywords=options.Keywords);
end
