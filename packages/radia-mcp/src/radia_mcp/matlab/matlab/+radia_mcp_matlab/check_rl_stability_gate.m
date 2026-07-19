function check_rl_stability_gate(seed)
r=radia_mcp_matlab.rl_stability_gate(seed);
disp(jsonencode(struct("tool","matlab_rl_stability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:RlStabilityGateFailed","RL stability gate failed."); end
end
