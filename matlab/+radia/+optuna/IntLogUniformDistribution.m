function distribution=IntLogUniformDistribution(low,high,step)
%INTLOGUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.
arguments
    low (1,1) double
    high (1,1) double
    step (1,1) double = 1
end
warning("radia:optuna:FutureWarning", ...
    "IntLogUniformDistribution is deprecated; use IntDistribution(Log=true).");
distribution=radia.optuna.IntDistribution(low,high,Log=true,Step=step);
distribution.name="IntLogUniformDistribution";
end
