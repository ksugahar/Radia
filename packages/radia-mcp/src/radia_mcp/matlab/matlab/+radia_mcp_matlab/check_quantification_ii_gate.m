function check_quantification_ii_gate(seed)
r=radia_mcp_matlab.quantification_ii_gate(seed); disp(jsonencode(struct("tool","matlab_quantification_ii_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:QuantificationIiGateFailed","Quantification II gate failed."); end
end
