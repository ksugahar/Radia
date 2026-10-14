function [A, info] = matrix_dump(vol_path, space, order, form, options)
%MATRIX_DUMP Assemble an NGSolve volume matrix and return sparse data.
%   Canonical snake_case MATLAB name matching the MEX/Python contract.

arguments
    vol_path (1,1) string
    space (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    form (1,1) string
    options.no_grads (1,1) logical = true
end

if nargout > 1
    [A, info] = radia.ngsolveMatrix(vol_path, space, order, form, ...
        NoGrads=options.no_grads);
else
    A = radia.ngsolveMatrix(vol_path, space, order, form, ...
        NoGrads=options.no_grads);
end
end
