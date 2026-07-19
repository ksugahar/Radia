function check_diffusion_inverse_gate(seed)
r=radia_mcp_matlab.diffusion_inverse_gate(seed);
disp(jsonencode(struct("tool","matlab_diffusion_inverse_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:DiffusionInverseGateFailed","Diffusion inverse gate failed."); end
end
