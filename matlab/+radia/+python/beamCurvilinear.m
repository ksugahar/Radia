function result = beamCurvilinear(functionName, positional, options)
%BEAMCURVILINEAR Call EarlyTimes curvilinear field-map operations.
%   This explicit batch boundary exposes Bishop/RMF loft-chain construction,
%   NGSolve-owned HCurl/HDiv projection, full transverse A sampling, and the
%   exact linear sample-to-polynomial-jet response. It is intended for model
%   initialization and artifact generation, not Simulink step-time use.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.beam_curvilinear",functionName,positional, ...
    Keywords=options.Keywords);
end
