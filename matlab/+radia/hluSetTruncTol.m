function hluSetTruncTol(tol)
%HLUSETTRUNCTOL Set HACApK H-LU low-rank recompression tolerance.
radia.internal.callMex('hlu.set_trunc_tol', double(tol));
end
