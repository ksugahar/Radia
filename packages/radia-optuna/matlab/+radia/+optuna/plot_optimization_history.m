function value=plot_optimization_history(study,varargin)
value=radia.optuna.internal.upstreamVisualization( ...
    "plot_optimization_history",study,varargin{:});
end
