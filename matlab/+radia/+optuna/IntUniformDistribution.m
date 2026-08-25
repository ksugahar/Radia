function distribution=IntUniformDistribution(low,high,step)
%INTUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.
arguments
    low (1,1) double
    high (1,1) double
    step (1,1) double = 1
end
warning("radia:optuna:FutureWarning", ...
    "IntUniformDistribution is deprecated; use IntDistribution.");
distribution=radia.optuna.IntDistribution(low,high,Step=step);
distribution.name="IntUniformDistribution";
end
