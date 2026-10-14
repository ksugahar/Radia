function check_sde_ito_gate(seed)
r=radia_mcp_matlab.sde_ito_gate(seed); disp(jsonencode(struct("tool","matlab_sde_ito_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:SdeItoGateFailed","SDE Ito gate failed."); end
end
