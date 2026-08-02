function result = streamFunction(functionName, positional, options)
%STREAMFUNCTION Call the Python stream-function design API explicitly.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.stream_function", ...
    functionName, positional, Keywords=options.Keywords);
end
