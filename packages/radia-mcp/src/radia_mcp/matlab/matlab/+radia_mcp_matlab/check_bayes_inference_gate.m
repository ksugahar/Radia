function check_bayes_inference_gate(seed)
r=radia_mcp_matlab.bayes_inference_gate(seed); disp(jsonencode(struct("tool","matlab_bayes_inference_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:BayesInferenceGateFailed","Bayes inference gate failed."); end
end
