function result = IC(matrix, options)
%IC Native real/complex incomplete Cholesky preconditioner.
arguments
    matrix (1,1) radia.ngsolve.Matrix
    options.Shift (1,1) double {mustBeFinite,mustBePositive} = 1.05
end
h = radia.internal.callMex('sparsesolv.ic', matrix.nativeHandle(), options.Shift);
result = radia.ngsolve.Matrix.fromNativeHandle(h);
end
