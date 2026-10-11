function result = ffagTaylorTopopt(functionName, positional, options)
%FFAGTAYLORTOPOPT Call the Python multi-momentum Taylor-map FFAG API.
%   This explicit batch boundary owns fixed-orbit R/T objective fusion and
%   HDiv-MMM topology optimization. It is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.ffag_taylor_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
