function result = sparsesolv(functionName, positional, options)
%SPARSESOLV Explicit Python batch boundary for the complete sparsesolv module.
%   Accepts Python NGSolve objects, never native MEX uint64 handles.
%   Use radia.sparsesolv.AMS/IC/COCR for the native matrix-handle path.
%   Remaining Python options, Update and solver-result access stay available
%   through the returned Python object. Not a Simulink per-step backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.sparsesolv_ngsolve", ...
    functionName, positional, Keywords=options.Keywords);
end
