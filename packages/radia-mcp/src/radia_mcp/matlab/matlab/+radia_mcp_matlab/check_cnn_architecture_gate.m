function check_cnn_architecture_gate(seed)
r=radia_mcp_matlab.cnn_architecture_gate(seed); disp(jsonencode(struct("tool","matlab_cnn_architecture_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:CnnArchitectureGateFailed","CNN architecture gate failed."); end
end
