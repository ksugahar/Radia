function [bRe, bIm] = bFromTrianglesComplex(vertices, currentRe, currentIm, obs, nThreads)
%BFROMTRIANGLESCOMPLEX Evaluate complex B from surface current triangles.

if nargin < 5
    nThreads = 0;
end
[bRe, bIm] = radia.internal.callMex( ...
    'biot_savart.b_triangles_complex', double(vertices), double(currentRe), ...
    double(currentIm), double(obs), double(nThreads));
end
