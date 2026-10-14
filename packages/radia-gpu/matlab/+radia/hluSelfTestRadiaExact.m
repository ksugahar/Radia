function err = hluSelfTestRadiaExact()
%HLUSELFTESTRADIAEXACT Run the fixed Radia-shape H-LU test.
err = radia.internal.callMex('hlu.self_test_radia_exact');
end
