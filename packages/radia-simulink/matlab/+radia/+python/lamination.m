function result = lamination(functionName, positional, options)
%LAMINATION Call the canonical Python lamination API explicitly.
%   The fallback is available for setup, validation, and batch computation;
%   it is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.lamination", functionName, positional, ...
    Keywords=options.Keywords);
end
