function result = ihAxisymCoupled(functionName, positional, options)
%IHAXISYMCOUPLED Call the canonical Python coupled axisymmetric IH solver.
%   Covers run_coupled (volumetric A_phi eddy currents with sigma(T),
%   mu_r(T) re-solved from the temperature, staggered with the enthalpy
%   heat integrator), AxisymEddyCurrent and EMMaterialTable.  The fallback
%   is available for setup, validation, and batch computation; it is not a
%   Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.ih_axisym_coupled", functionName, positional, ...
    Keywords=options.Keywords);
end