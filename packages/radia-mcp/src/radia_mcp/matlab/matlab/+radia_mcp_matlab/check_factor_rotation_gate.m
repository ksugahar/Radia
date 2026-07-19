function check_factor_rotation_gate(seed)
r=radia_mcp_matlab.factor_rotation_gate(seed); disp(jsonencode(struct("tool","matlab_factor_rotation_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:FactorRotationGateFailed","Factor rotation gate failed."); end
end
