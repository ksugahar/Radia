function err = hluSelfTestDepth3Asymmetric(nbTiny)
%HLUSELFTESTDEPTH3ASYMMETRIC Run the asymmetric depth-three test.
if nargin < 1, nbTiny = 3; end
err = radia.internal.callMex('hlu.self_test_depth3_asymmetric', double(nbTiny));
end
