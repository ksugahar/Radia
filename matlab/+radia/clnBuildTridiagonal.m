function T = clnBuildTridiagonal(diagValues)
%CLNBUILDTRIDIAGONAL Build the CLN tridiagonal inductance matrix.
T = radia.internal.callMex('cln.build_tridiagonal', double(diagValues));
end
