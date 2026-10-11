function result = axisymMeasure(functionName, positional, options)
% Explicit batch entry for radia.axisym_measure; not a per-step backend.
% An axisymmetric mesh carries two integration measures that look identical
% in code -- the meridian dr dz and the revolved 2 pi r dr dz -- so the
% Python side names them rather than leaving the choice to memory, and this
% entry keeps that naming intact for MATLAB callers.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.axisym_measure", ...
    functionName, positional, Keywords=options.Keywords);
end
