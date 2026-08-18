function result = beamYoshidaSplit(functionName, positional, options)
%BEAMYOSHIDASPLIT Call Yoshida split-operator tracking on a chain field.
%   Explicit batch boundary over radia.beam_yoshida_split: the original
%   EarlyTimes symmetric split-operator integrator (Ishi 2016 lineage)
%   revived on the CanonicalHCurl chain polynomials -- exact symplectic
%   factors at any amplitude, order 2 or 4 (Yoshida w-coefficients),
%   paraxial Hamiltonian. Intended for model initialization and
%   verification tracking, not Simulink step-time use.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.beam_yoshida_split",functionName,positional, ...
    Keywords=options.Keywords);
end
