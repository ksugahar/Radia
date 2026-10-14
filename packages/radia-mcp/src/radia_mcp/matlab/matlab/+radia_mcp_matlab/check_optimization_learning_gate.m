function check_optimization_learning_gate(seed)
r=radia_mcp_matlab.optimization_learning_gate(seed);
disp(jsonencode(struct("tool","matlab_optimization_learning_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:OptimizationLearningGateFailed","Optimization learning gate failed."); end
end
