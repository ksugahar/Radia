function [hRe, hIm] = hFromSegmentsComplex(segments, obs, currentRe, currentIm, nThreads)
%HFROMSEGMENTSCOMPLEX Evaluate complex H from finite current segments.

if nargin < 5
    nThreads = 0;
end
[hRe, hIm] = radia.internal.callMex( ...
    'biot_savart.h_segments_complex', double(segments), double(obs), ...
    double(currentRe), double(currentIm), double(nThreads));
end
