function result = machines(moduleName, functionName, positional, options)
%MACHINES Call a Python-only motor, MagLev, PCB, or stream-function operation.
arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
allowed = any(startsWith(moduleName, ["motor_","maglev.","streamfunction"])) || ...
    ismember(moduleName, ["magnet","pcb_design","stream_function"]);
if ~allowed
    error("radia:python:Module", "Unsupported machine/application module: %s", moduleName);
end
result = radia.internal.callPython("radia." + moduleName, functionName, ...
    positional, Keywords=options.Keywords);
end
