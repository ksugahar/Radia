function check_ml_cae_pipeline(stage, seed)
%CHECK_ML_CAE_PIPELINE MCP-facing JSON wrapper for the CAE/ML pipeline.
report = radia_mcp_matlab.ml_cae_pipeline(string(stage), seed);
disp(jsonencode(struct("tool", "matlab_ml_cae_pipeline", ...
    "ok", report.ok, "result", report)));
if ~report.ok
    error("radia_mcp_matlab:MlCaeGateFailed", ...
        "ML/CAE stage '%s' did not pass its gate.", string(stage));
end
end
