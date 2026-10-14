function check_linear_inverse_gate(seed)
r=radia_mcp_matlab.linear_inverse_gate(seed); disp(jsonencode(struct("tool","matlab_linear_inverse_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:LinearInverseGateFailed","Linear inverse gate failed."); end
end
