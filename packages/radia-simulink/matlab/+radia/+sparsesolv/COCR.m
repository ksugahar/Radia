function result = COCR(matrix, preconditioner, options)
%COCR Native inverse operator for symmetric (not Hermitian) systems.
%   Apply with result.matvec(rhs). The operator retains its matrix and
%   preconditioner even after their MATLAB owners are deleted.
arguments
    matrix (1,1) radia.ngsolve.Matrix
    preconditioner (1,1) radia.ngsolve.Matrix
    options.Tolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-8
    options.MaxSteps (1,1) double {mustBeInteger,mustBePositive} = 500
end
h = radia.internal.callMex('sparsesolv.cocr', matrix.nativeHandle(), ...
    preconditioner.nativeHandle(), options.Tolerance, options.MaxSteps);
result = radia.ngsolve.Matrix.fromNativeHandle(h);
end
