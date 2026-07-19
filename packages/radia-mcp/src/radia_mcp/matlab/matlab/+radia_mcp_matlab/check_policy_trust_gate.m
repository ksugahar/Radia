function check_policy_trust_gate(seed)
r=radia_mcp_matlab.policy_trust_gate(seed); disp(jsonencode(struct("tool","matlab_policy_trust_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:PolicyTrustGateFailed","Policy trust gate failed."); end
end
