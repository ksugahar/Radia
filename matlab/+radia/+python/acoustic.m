function result = acoustic(moduleName, functionName, positional, options)
%ACOUSTIC Call a Python-only acoustic CQ or FSI operation.
arguments
    moduleName (1,1) string {mustBeMember(moduleName,["cq","fsi"])}
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.acoustics." + moduleName, ...
    functionName, positional, Keywords=options.Keywords);
end
