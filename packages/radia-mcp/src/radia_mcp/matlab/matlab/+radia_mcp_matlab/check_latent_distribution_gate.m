function check_latent_distribution_gate(seed)
r=radia_mcp_matlab.latent_distribution_gate(seed); disp(jsonencode(struct("tool","matlab_latent_distribution_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:LatentDistributionGateFailed","Latent distribution gate failed."); end
end
