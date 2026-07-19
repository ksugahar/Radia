function check_neural_ode_gate(seed)
r=radia_mcp_matlab.neural_ode_gate(seed); disp(jsonencode(struct("tool","matlab_neural_ode_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:NeuralOdeGateFailed","Neural ODE gate failed."); end
end
