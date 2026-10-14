function info = space_info(vol_path, order, options)
%SPACE_INFO Return NGSolve finite-element space dimensions.
%   Canonical snake_case MATLAB name matching the MEX/Python contract.

arguments
    vol_path (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    options.no_grads (1,1) logical = true
end

info = radia.spaceInfo(vol_path, order, NoGrads=options.no_grads);
end
