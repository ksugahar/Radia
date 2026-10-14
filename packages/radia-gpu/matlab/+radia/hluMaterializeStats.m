function stats = hluMaterializeStats()
%HLUMATERIALIZESTATS Return H-LU materialization fallback counters.
stats = radia.internal.callMex('hlu.materialize_stats');
end
