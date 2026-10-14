function err = hluSelfTestRadiaExactWithMatrix(A, b)
%HLUSELFTESTRADIAEXACTWITHMATRIX Run the 162-by-162 Radia-shape test.
err = radia.internal.callMex('hlu.self_test_radia_exact_with_matrix', double(A), double(b));
end
