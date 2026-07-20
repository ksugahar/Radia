function err = hluSelfTestRk(nPerBlock, rkRank)
%HLUSELFTESTRK Run the rank-aware HACApK H-LU self-test.
if nargin < 1, nPerBlock = 100; end
if nargin < 2, rkRank = 5; end
err = radia.internal.callMex('hlu.self_test_rk', double(nPerBlock), double(rkRank));
end
