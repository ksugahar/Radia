function result = p1Linear(functionName, positional, options)
% Explicit batch entry for radia.p1_linear; not a per-step backend.
% Linear total-A magnetostatics with one lowest-order beta-zero AMS-PCG solve.
% Mesh, GridFunction and CoefficientFunction objects remain on the Python side.
% Call outside Python TaskManager: this solver owns assembly/solve regions.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.p1_linear", ...
    functionName, positional, Keywords=options.Keywords);
end
