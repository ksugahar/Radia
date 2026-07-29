function result = isochronousTopopt(functionName, positional, options)
%ISOCHRONOUSTOPOPT Call the Python HDiv-MMM density-topology API.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.isochronous_topopt", ...
    functionName, positional, Keywords=options.Keywords);
end
