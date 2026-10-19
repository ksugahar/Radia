function [aRe, aIm] = aFromSegmentsComplex(segments, obs, currentRe, currentIm, nThreads)
%AFROMSEGMENTSCOMPLEX Evaluate complex A from finite current segments.

if nargin < 5
    nThreads = 0;
end
[aRe, aIm] = radia.internal.callMex( ...
    'biot_savart.a_segments_complex', double(segments), double(obs), ...
    double(currentRe), double(currentIm), double(nThreads));
end
