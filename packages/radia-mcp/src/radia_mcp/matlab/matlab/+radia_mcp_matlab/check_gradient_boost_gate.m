function check_gradient_boost_gate(seed)
r=radia_mcp_matlab.gradient_boost_gate(seed); disp(jsonencode(struct("tool","matlab_gradient_boost_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:GradientBoostGateFailed","Gradient boost gate failed."); end
end
