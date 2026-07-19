function result = clnLanczos(K, N, nIter, tol)
%CLNLANCZOS Reduce two PEEC matrices to the CLN I form.
arguments
    K double
    N double
    nIter (1,1) double = -1
    tol (1,1) double = 1e-30
end
result = radia.internal.callMex('cln.lanczos', double(K), double(N), nIter, tol);
end
