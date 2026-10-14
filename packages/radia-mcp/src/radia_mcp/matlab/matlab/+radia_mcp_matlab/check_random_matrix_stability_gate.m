function check_random_matrix_stability_gate(seed)
r=radia_mcp_matlab.random_matrix_stability_gate(seed); disp(jsonencode(struct("tool","matlab_random_matrix_stability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:RandomMatrixStabilityGateFailed","Random matrix stability gate failed."); end
end
