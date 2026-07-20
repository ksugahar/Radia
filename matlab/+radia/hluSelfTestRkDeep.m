function err = hluSelfTestRkDeep(nPerBlock, rkRank)
%HLUSELFTESTRKDEEP Run the deep rank-aware HACApK H-LU self-test.
if nargin < 1, nPerBlock = 100; end
if nargin < 2, rkRank = 5; end
err = radia.internal.callMex('hlu.self_test_rk_deep', double(nPerBlock), double(rkRank));
end
