function check_dueling_cae_action_gate(seed)
r=radia_mcp_matlab.dueling_cae_action_gate(seed);
disp(jsonencode(struct("tool","matlab_dueling_cae_action_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:DuelingCaeActionGateFailed","Dueling CAE action gate failed."); end
end
