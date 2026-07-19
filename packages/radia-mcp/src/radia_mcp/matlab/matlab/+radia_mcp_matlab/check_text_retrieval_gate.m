function check_text_retrieval_gate(seed)
r=radia_mcp_matlab.text_retrieval_gate(seed); disp(jsonencode(struct("tool","matlab_text_retrieval_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TextRetrievalGateFailed","Text retrieval gate failed."); end
end
