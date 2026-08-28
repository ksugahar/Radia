function distribution=DiscreteUniformDistribution(low,high,q)
%DISCRETEUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.
warning("radia:optuna:FutureWarning", ...
    "DiscreteUniformDistribution is deprecated; use FloatDistribution(Step=q).");
distribution=radia.optuna.FloatDistribution(low,high,Step=q);
distribution.name="DiscreteUniformDistribution";
end
