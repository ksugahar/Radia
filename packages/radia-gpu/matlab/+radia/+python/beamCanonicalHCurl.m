function result = beamCanonicalHCurl(functionName, positional, options)
%BEAMCANONICALHCURL Call CanonicalHCurl vacuum-chain operations.
%   Explicit batch boundary over radia.beam_canonical_hcurl: the
%   tracking-specialized vacuum HCurl subspace (dimension law
%   p_x*(p_s+1), gauge-rigid, curved-chart least-defect basis), its
%   full-volume frame-B projection, and the per-segment covariant
%   transverse coefficient arrays consumed by the Lie kernel. Intended
%   for model initialization and artifact generation, not Simulink
%   step-time use.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.beam_canonical_hcurl",functionName,positional, ...
    Keywords=options.Keywords);
end
