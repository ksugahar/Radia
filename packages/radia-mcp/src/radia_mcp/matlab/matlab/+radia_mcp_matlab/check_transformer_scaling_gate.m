function check_transformer_scaling_gate(seed)
r=radia_mcp_matlab.transformer_scaling_gate(seed); disp(jsonencode(struct("tool","matlab_transformer_scaling_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TransformerScalingGateFailed","Transformer scaling gate failed."); end
end
