function err = hluSelfTestRadiaExactDiag(diagBoost)
%HLUSELFTESTRADIAEXACTDIAG Run the adjustable diagonal-dominance test.
if nargin < 1, diagBoost = 2; end
err = radia.internal.callMex('hlu.self_test_radia_exact_diag', double(diagBoost));
end
