function check_activation_probability_gate(seed)
r=radia_mcp_matlab.activation_probability_gate(seed); disp(jsonencode(struct("tool","matlab_activation_probability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:ActivationProbabilityGateFailed","Activation probability gate failed."); end
end
