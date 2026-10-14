function check_gan_design_gate(seed)
r=radia_mcp_matlab.gan_design_gate(seed);
disp(jsonencode(struct("tool","matlab_gan_design_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:GanDesignGateFailed","GAN design gate failed."); end
end
