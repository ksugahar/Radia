function result = acceleratorMagnetTopopt(functionName, positional, options)
%ACCELERATORMAGNETTOPOPT Call the Python design-orbit HDiv-MMM API.
%   This explicit batch boundary is intended for setup, optimization runs,
%   and artifact generation. It is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_magnet_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
