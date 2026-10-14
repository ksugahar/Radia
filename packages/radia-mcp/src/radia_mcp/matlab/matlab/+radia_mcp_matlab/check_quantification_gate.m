function check_quantification_gate(seed)
r=radia_mcp_matlab.quantification_gate(seed); disp(jsonencode(struct("tool","matlab_quantification_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:QuantificationGateFailed","Quantification gate failed."); end
end
