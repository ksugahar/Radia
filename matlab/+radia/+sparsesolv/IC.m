function result = IC(matrix, options)
%IC Native real/complex incomplete Cholesky preconditioner.
%   Uses transpose IC for real SPD or complex-symmetric matrices. For
%   Hermitian positive-definite systems use ICCG with Conjugate=true.
%   Shift is the multiplicative IC diagonal shift alpha (alpha*a_ii on rows
%   with Re(a_ii) > 0) and must be >= 1: the native factorization rejects a
%   smaller value with radia:mex:Exception (this wrapper only checks that it
%   is positive). The shift is fixed, with no automatic search, and a zero or
%   non-finite pivot also raises radia:mex:Exception.
arguments
    matrix (1,1) radia.ngsolve.Matrix
    options.Shift (1,1) double {mustBeFinite,mustBePositive} = 1.05
end
h = radia.internal.callMex('sparsesolv.ic', matrix.nativeHandle(), options.Shift);
result = radia.ngsolve.Matrix.fromNativeHandle(h);
end
