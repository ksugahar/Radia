function [singleLayer, doubleLayer] = assembleSldlGalerkin( ...
        vertices, triangles, p2Nodes, regularDegree, singularOrder, nThreads)
%ASSEMBLESLDLGALERKIN Assemble P1 Galerkin Laplace SL/DL matrices.

if nargin < 4
    regularDegree = 11;
end
if nargin < 5
    singularOrder = 8;
end
if nargin < 6
    nThreads = 0;
end
[singleLayer, doubleLayer] = radia.internal.callMex( ...
    'bem.assemble_sldl', double(vertices), int64(triangles), ...
    double(p2Nodes), double(regularDegree), double(singularOrder), ...
    double(nThreads));
end
