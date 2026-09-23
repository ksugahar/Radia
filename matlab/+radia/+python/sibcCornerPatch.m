function result = sibcCornerPatch(functionName, positional, options)
% Explicit batch entry for radia.sibc_corner_patch; not a per-step backend.
% The matched patch that repairs a surface impedance near tips and corners,
% and the whole-section control it is scored against. Batch only: a patch
% owns an NGSolve mesh and GridFunction, which stay on the Python side.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.sibc_corner_patch", ...
    functionName, positional, Keywords=options.Keywords);
end