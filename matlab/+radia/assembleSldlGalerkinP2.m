function [singleLayer, doubleLayer] = assembleSldlGalerkinP2( ...
        vertices, triangles, p2Nodes, dofsPerTriangle, nDof, ...
        regularDegree, singularOrder, nThreads)
%ASSEMBLESLDLGALERKINP2 Assemble P2 Lagrange Galerkin Laplace SL/DL matrices.

if nargin < 6
    regularDegree = 11;
end
if nargin < 7
    singularOrder = 8;
end
if nargin < 8
    nThreads = 0;
end
[singleLayer, doubleLayer] = radia.internal.callMex( ...
    'bem.assemble_sldl_p2', double(vertices), int64(triangles), ...
    double(p2Nodes), int64(dofsPerTriangle), double(nDof), ...
    double(regularDegree), double(singularOrder), double(nThreads));
end
