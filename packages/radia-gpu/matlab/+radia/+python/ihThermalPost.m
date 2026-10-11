function result = ihThermalPost(functionName, positional, options)
%IHTHERMALPOST Call the canonical Python IH temperature post-processing API.
%   Covers thermal_exposure (volume above temperature limits, hottest-ring
%   spread) and case_depth (threshold depth below the heated surface).
%   The fallback is available for setup, validation, and batch
%   computation; it is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.ih_thermal_post", functionName, positional, ...
    Keywords=options.Keywords);
end