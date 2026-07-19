function check_time_series_gate(seed)
r=radia_mcp_matlab.time_series_gate(seed); disp(jsonencode(struct("tool","matlab_time_series_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TimeSeriesGateFailed","Time series gate failed."); end
end
