function check_tabular_transform_gate(seed)
r=radia_mcp_matlab.tabular_transform_gate(seed); disp(jsonencode(struct("tool","matlab_tabular_transform_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TabularTransformGateFailed","Tabular transform gate failed."); end
end
