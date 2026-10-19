function err = hluSelfTestMixedSiblingViaConversion(nbSmall)
%HLUSELFTESTMIXEDSIBLINGVIACONVERSION Run the HACApK conversion-path test.
if nargin < 1, nbSmall = 5; end
err = radia.internal.callMex('hlu.self_test_mixed_sibling_via_conversion', double(nbSmall));
end
