function result = clnAcaCompress(P, eps, kMax)
%CLNACACOMPRESS Compute a low-rank ACA/SVD compression of a star matrix.
arguments
    P double
    eps (1,1) double = 1e-4
    kMax (1,1) double = -1
end
result = radia.internal.callMex('cln.aca_compress', double(P), eps, kMax);
end
