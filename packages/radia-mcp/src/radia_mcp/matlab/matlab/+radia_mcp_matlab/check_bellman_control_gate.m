function check_bellman_control_gate(seed)
r=radia_mcp_matlab.bellman_control_gate(seed); disp(jsonencode(struct("tool","matlab_bellman_control_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:BellmanControlGateFailed","Bellman control gate failed."); end
end
