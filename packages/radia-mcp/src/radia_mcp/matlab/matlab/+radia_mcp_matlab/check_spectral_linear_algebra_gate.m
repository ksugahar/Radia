function check_spectral_linear_algebra_gate(seed)
r=radia_mcp_matlab.spectral_linear_algebra_gate(seed); disp(jsonencode(struct("tool","matlab_spectral_linear_algebra_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:SpectralLinearAlgebraGateFailed","Spectral linear algebra gate failed."); end
end
