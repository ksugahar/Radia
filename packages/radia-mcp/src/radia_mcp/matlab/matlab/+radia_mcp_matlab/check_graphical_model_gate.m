function check_graphical_model_gate(seed)
r=radia_mcp_matlab.graphical_model_gate(seed); disp(jsonencode(struct("tool","matlab_graphical_model_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:GraphicalModelGateFailed","Graphical model gate failed."); end
end
