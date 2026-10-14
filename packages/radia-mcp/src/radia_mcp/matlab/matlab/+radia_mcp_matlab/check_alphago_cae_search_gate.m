function check_alphago_cae_search_gate(seed)
r=radia_mcp_matlab.alphago_cae_search_gate(seed);
disp(jsonencode(struct("tool","matlab_alphago_cae_search_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:AlphaGoCaeSearchGateFailed","AlphaGo CAE search gate failed."); end
end
