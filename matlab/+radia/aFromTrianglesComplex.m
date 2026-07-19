function [aRe, aIm] = aFromTrianglesComplex(vertices, currentRe, currentIm, obs, nThreads)
%AFROMTRIANGLESCOMPLEX Evaluate complex A from surface current triangles.

if nargin < 5
    nThreads = 0;
end
[aRe, aIm] = radia.internal.callMex( ...
    'biot_savart.a_triangles_complex', double(vertices), double(currentRe), ...
    double(currentIm), double(obs), double(nThreads));
end
