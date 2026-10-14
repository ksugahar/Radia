function check_variational_geometry_gate(seed)
r=radia_mcp_matlab.variational_geometry_gate(seed); disp(jsonencode(struct("tool","matlab_variational_geometry_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:VariationalGeometryGateFailed","Variational geometry gate failed."); end
end
