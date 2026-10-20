function result = acceleratorFieldValidation(functionName, positional, options)
%ACCELERATORFIELDVALIDATION Call the Python observation-tube API explicitly.
%   This batch boundary compares gauge-invariant accelerator field samples.
%   It is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_field_validation", functionName, positional, ...
    Keywords=options.Keywords);
end
