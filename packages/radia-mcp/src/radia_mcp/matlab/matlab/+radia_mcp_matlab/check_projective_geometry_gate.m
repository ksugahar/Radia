function check_projective_geometry_gate(seed)
r=radia_mcp_matlab.projective_geometry_gate(seed); disp(jsonencode(struct("tool","matlab_projective_geometry_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:ProjectiveGeometryGateFailed","Projective geometry gate failed."); end
end
