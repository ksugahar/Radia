function distribution=FloatDistribution(low,high,options)
%FLOATDISTRIBUTION Public Optuna-style float distribution descriptor.
arguments
    low (1,1) double
    high (1,1) double
    options.Log (1,1) logical = false
    options.Step (1,1) double = NaN
end
[effectiveHigh,adjusted]=radia.optuna.internal.UpstreamNumerics. ...
    adjustDiscreteUniformHigh(low,high,options.Step);
if adjusted
    warning("radia:optuna:DistributionAdjusted", ...
        "FloatDistribution high was adjusted from %g to %g to align with Step=%g.", ...
        high,effectiveHigh,options.Step);
end
distribution=radia.optuna.internal.DistributionCodec.float( ...
    low,effectiveHigh,options.Log,options.Step);
distribution.name="FloatDistribution";
end
