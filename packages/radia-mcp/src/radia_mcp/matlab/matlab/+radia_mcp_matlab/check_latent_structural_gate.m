function check_latent_structural_gate(seed)
r=radia_mcp_matlab.latent_structural_gate(seed); disp(jsonencode(struct("tool","matlab_latent_structural_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:LatentStructuralGateFailed","Latent structural gate failed."); end
end
