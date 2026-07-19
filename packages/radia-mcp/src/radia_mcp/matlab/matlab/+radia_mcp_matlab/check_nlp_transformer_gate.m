function check_nlp_transformer_gate(seed)
r=radia_mcp_matlab.nlp_transformer_gate(seed); disp(jsonencode(struct("tool","matlab_nlp_transformer_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:NlpTransformerGateFailed","NLP transformer gate failed."); end
end
