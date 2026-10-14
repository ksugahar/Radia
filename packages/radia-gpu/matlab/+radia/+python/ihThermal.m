function result = ihThermal(functionName, positional, options)
%IHTHERMAL Call the canonical Python induction-heating thermal API explicitly.
%   Covers field-artifact sidecars (load_field, write_field_sidecar), the
%   verified EM-to-thermal heat-source transfer and the exact
%   circumferential average.  The fallback is available for setup,
%   validation, and batch computation; it is not a Simulink step-time
%   backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.ih_thermal", functionName, positional, ...
    Keywords=options.Keywords);
end