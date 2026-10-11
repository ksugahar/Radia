function tol = hluGetTruncTol()
%HLUGETTRUNCTOL Return HACApK H-LU low-rank recompression tolerance.
tol = radia.internal.callMex('hlu.get_trunc_tol');
end
