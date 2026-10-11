function result = meshedCurrent(functionName, positional, options)
% Explicit batch entry for radia.meshed_current; not a per-step backend.
% Impressed DC current of a closed meshed conductor: the mixed RT0/P0 current
% and the A-phi scalar-potential current with one cut. Batch only: the
% current is an NGSolve GridFunction on the conductor mesh, which stays on the
% Python side. Python ecosystem objects remain in result.value; native MEX
% handles cannot be passed to this interface.
% Parallel assembly uses a caller-owned TaskManager.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.meshed_current", ...
    functionName, positional, Keywords=options.Keywords);
end
