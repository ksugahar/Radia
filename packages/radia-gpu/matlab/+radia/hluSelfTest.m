function err = hluSelfTest(depth, nPerBlock)
%HLUSELFTEST Run the dense-leaf HACApK H-LU self-test.
if nargin < 1, depth = 1; end
if nargin < 2, nPerBlock = 100; end
err = radia.internal.callMex('hlu.self_test', double(depth), double(nPerBlock));
end
