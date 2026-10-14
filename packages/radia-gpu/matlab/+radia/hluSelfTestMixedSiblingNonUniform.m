function err = hluSelfTestMixedSiblingNonUniform(n1, n2, m1, m3)
%HLUSELFTESTMIXEDSIBLINGNONUNIFORM Run the nonuniform mixed-sibling test.
if nargin < 1, n1 = 5; end
if nargin < 2, n2 = 7; end
if nargin < 3, m1 = 2; end
if nargin < 4, m3 = 3; end
err = radia.internal.callMex('hlu.self_test_mixed_sibling_nonuniform', ...
    double(n1), double(n2), double(m1), double(m3));
end
