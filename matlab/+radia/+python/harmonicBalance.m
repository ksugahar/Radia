function result = harmonicBalance(functionName, positional, options)
%HARMONICBALANCE Call the Python odd-harmonic balance batch API.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.harmonic_balance", ...
    functionName, positional, Keywords=options.Keywords);
end
