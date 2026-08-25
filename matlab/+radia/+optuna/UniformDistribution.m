function distribution=UniformDistribution(low,high)
%UNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.
warning("radia:optuna:FutureWarning", ...
    "UniformDistribution is deprecated; use FloatDistribution.");
distribution=radia.optuna.FloatDistribution(low,high);
distribution.name="UniformDistribution";
end
