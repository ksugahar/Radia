function value=is_available(options)
arguments
    options.Backend (1,1) string = "plotly"
end
value=radia.optuna.internal.upstreamVisualization( ...
    "is_available",[],"Backend",options.Backend);
end
