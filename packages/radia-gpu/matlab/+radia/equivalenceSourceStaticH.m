function H = equivalenceSourceStaticH(centroids, normals, areas, Hsurf, obs, nThreads)
%EQUIVALENCESOURCESTATICH Reconstruct exterior magnetostatic H from an equivalent surface.
if nargin < 6, nThreads = 0; end
H = radia.internal.callMex('equivalence.static_h', double(centroids), ...
    double(normals), double(areas), double(Hsurf), double(obs), double(nThreads));
end
