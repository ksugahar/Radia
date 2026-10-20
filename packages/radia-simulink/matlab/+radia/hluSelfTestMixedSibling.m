function err = hluSelfTestMixedSibling(nbSmall)
%HLUSELFTESTMIXEDSIBLING Run the mixed-sibling H-LU test.
if nargin < 1, nbSmall = 5; end
err = radia.internal.callMex('hlu.self_test_mixed_sibling', double(nbSmall));
end
