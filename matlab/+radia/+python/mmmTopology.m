function result = mmmTopology(functionName, positional, options)
%MMMTOPOLOGY Call the two-stage HDiv-MMM topology API explicitly.
%   This named batch fallback preserves the canonical binary-Lego then
%   GetTrafo ordering. It is not permitted inside a Simulink time step.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.mmm_topology", ...
    functionName, positional, Keywords=options.Keywords);
end
