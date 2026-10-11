function result = xsuiteBridge(functionName, positional, options)
%XSUITEBRIDGE Call the Python Radia-to-Xsuite batch bridge.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.xsuite_bridge", ...
    functionName, positional, Keywords=options.Keywords);
end
