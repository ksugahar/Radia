function check_pca_covariance_gate(seed)
r=radia_mcp_matlab.pca_covariance_gate(seed); disp(jsonencode(struct("tool","matlab_pca_covariance_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:PcaCovarianceGateFailed","PCA covariance gate failed."); end
end
