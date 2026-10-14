function check_integrity_encoding_gate(seed)
r=radia_mcp_matlab.integrity_encoding_gate(seed); disp(jsonencode(struct("tool","matlab_integrity_encoding_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:IntegrityEncodingGateFailed","Integrity encoding gate failed."); end
end
