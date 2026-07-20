function solverConfig(options)
%SOLVERCONFIG Set native Radia solver options from a scalar struct.

if nargin == 0
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error("radia:solverConfig:Options", "options must be a scalar struct.");
end
radia.internal.callMex('radia.SolverConfig', options);
end
