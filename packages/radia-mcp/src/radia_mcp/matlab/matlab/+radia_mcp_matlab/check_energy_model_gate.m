function check_energy_model_gate(seed)
r=radia_mcp_matlab.energy_model_gate(seed);
disp(jsonencode(struct("tool","matlab_energy_model_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:EnergyModelGateFailed","Energy model gate failed."); end
end
