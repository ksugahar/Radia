function config = getSolverConfig()
%GETSOLVERCONFIG Return the native Radia solver configuration.

config = radia.internal.callMex('radia.GetSolverConfig');
end
