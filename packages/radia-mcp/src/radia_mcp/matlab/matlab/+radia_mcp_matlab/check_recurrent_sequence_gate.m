function check_recurrent_sequence_gate(seed)
r=radia_mcp_matlab.recurrent_sequence_gate(seed); disp(jsonencode(struct("tool","matlab_recurrent_sequence_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:RecurrentSequenceGateFailed","Recurrent sequence gate failed."); end
end
