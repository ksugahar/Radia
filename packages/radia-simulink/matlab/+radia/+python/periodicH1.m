function result = periodicH1(functionName, positional, options)
% Explicit batch entry for radia.periodic_h1; not a per-step backend.
% Builds a face-consistent H1 space across one explicit periodic interface.
% The returned space is a Python ecosystem object and stays in result.value:
% it is an argument to later Python calls, not a value MATLAB can compute on,
% and it must be rebuilt after loading a mesh rather than serialized.
% Native MEX handles cannot be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.periodic_h1", ...
    functionName, positional, Keywords=options.Keywords);
end
