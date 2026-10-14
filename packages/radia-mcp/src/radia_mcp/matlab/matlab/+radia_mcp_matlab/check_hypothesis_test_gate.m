function check_hypothesis_test_gate(seed)
r=radia_mcp_matlab.hypothesis_test_gate(seed); disp(jsonencode(struct("tool","matlab_hypothesis_test_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:HypothesisTestGateFailed","Hypothesis test gate failed."); end
end
