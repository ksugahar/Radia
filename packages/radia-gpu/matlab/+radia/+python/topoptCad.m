function result = topoptCad(functionName, positional, options)
%TOPOPTCAD Call the Python topology-optimization CAD bridge.
%   Nodal level sets from per-element densities, Exodus level-set
%   hand-off for the Cubit `create tri iso` route, and marching-cubes
%   iso-surface STL (radia.topopt_cad).  The mesh-side half (STL ->
%   gated hex/tet .vol) is the cubit MCP tool `cubit_stl_to_vol`.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.topopt_cad", ...
    functionName, positional, Keywords=options.Keywords);
end
