function check_a3c_async_gate(seed)
r=radia_mcp_matlab.a3c_async_gate(seed);
disp(jsonencode(struct("tool","matlab_a3c_async_gate","ok",r.ok,"result",r)));
end
