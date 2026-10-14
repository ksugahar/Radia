function check_nice_flow_gate(seed)
r=radia_mcp_matlab.nice_flow_gate(seed);
disp(jsonencode(struct("tool","matlab_nice_flow_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:NiceFlowGateFailed","NICE flow gate failed."); end
end
