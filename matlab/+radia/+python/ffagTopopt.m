function result = ffagTopopt(functionName, positional, options)
%FFAGTOPOPT Call the Python FFAG transfer-matrix topology API.
%   This explicit batch boundary is intended for setup, optimization runs,
%   and artifact generation. It is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.ffag_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
