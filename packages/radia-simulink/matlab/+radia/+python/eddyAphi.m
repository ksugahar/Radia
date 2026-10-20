function result = eddyAphi(functionName, positional, options)
% Explicit batch entry for radia.eddy_aphi; not a per-step backend.
% Terminal-driven three-dimensional A-V eddy solves. Batch only: the solve
% owns an NGSolve mesh and GridFunction, which stay on the Python side.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.eddy_aphi", ...
    functionName, positional, Keywords=options.Keywords);
end