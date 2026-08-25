function distribution=LogUniformDistribution(low,high)
%LOGUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.
warning("radia:optuna:FutureWarning", ...
    "LogUniformDistribution is deprecated; use FloatDistribution(Log=true).");
distribution=radia.optuna.FloatDistribution(low,high,Log=true);
distribution.name="LogUniformDistribution";
end
