function info = spaceInfo(volPath, order, options)
%SPACEINFO Construct NGSolve HCurl/HDiv spaces and return their dimensions.

arguments
    volPath (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    options.NoGrads (1,1) logical = true
end

info = radia.internal.callMex( ...
    'ngsolve.space_info', char(volPath), order, options.NoGrads);
end
