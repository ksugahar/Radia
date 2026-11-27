function result = topoptCad(functionName, positional, options)
%TOPOPTCAD Call the Python topology-optimization CAD bridge.
%   Nodal level sets from per-element densities, Exodus level-set
%   hand-off for the Cubit `create tri iso` route, and marching-cubes
%   iso-surface STL (radia.topopt_cad).  Meshing the STL into a checked
%   hex/tet .vol is done by the caller outside this bridge.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.topopt_cad", ...
    functionName, positional, Keywords=options.Keywords);
end
