function check_tabular_lookup_gate(seed)
r=radia_mcp_matlab.tabular_lookup_gate(seed); disp(jsonencode(struct("tool","matlab_tabular_lookup_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TabularLookupGateFailed","Tabular lookup gate failed."); end
end
