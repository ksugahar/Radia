function check_vae_field_gate(seed)
r=radia_mcp_matlab.vae_field_gate(seed);
disp(jsonencode(struct("tool","matlab_vae_field_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:VaeFieldGateFailed","VAE field gate failed."); end
end
