function result = urn(moduleName, functionName, positional, options)
%URN Call a Python Universal Relaxation Network operation.
arguments
    moduleName (1,1) string {mustBeMember(moduleName, ...
        ["relaxation_basis_library","universal_relaxation_network"])}
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.urn." + moduleName, ...
    functionName, positional, Keywords=options.Keywords);
end
