function result = topologyOptimization(functionName, positional, options)
%TOPOLOGYOPTIMIZATION Call the Python HDiv-MMM topology API explicitly.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.topology_optimization", ...
    functionName, positional, Keywords=options.Keywords);
end
