function check_quantum_linear_gate(seed)
r=radia_mcp_matlab.quantum_linear_gate(seed); disp(jsonencode(struct("tool","matlab_quantum_linear_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:QuantumLinearGateFailed","Quantum linear gate failed."); end
end
