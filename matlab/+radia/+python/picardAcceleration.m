function result = picardAcceleration(functionName, positional, options)
%PICARDACCELERATION Explicit batch access to the shared Picard implementation.
% ConstrainedAndersonAccelerator is retained as a Python object in result.value;
% step, reset and stats operate on that same object. Not a Simulink step backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.picard_acceleration", ...
    functionName, positional, Keywords=options.Keywords);
end
