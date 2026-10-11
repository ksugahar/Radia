function err = hluSelfTestAddmulRkRk(m, n, inner, kA, kB, kC)
%HLUSELFTESTADDMULRKRK Run the rank-times-rank update self-test.
if nargin < 1, m = 64; end
if nargin < 2, n = 64; end
if nargin < 3, inner = 64; end
if nargin < 4, kA = 5; end
if nargin < 5, kB = 5; end
if nargin < 6, kC = 5; end
err = radia.internal.callMex('hlu.self_test_addmul_rkrk', ...
    double(m), double(n), double(inner), double(kA), double(kB), double(kC));
end
