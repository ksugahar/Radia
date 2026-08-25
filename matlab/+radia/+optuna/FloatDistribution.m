function distribution=FloatDistribution(low,high,options)
%FLOATDISTRIBUTION Public Optuna-style float distribution descriptor.
arguments
    low (1,1) double
    high (1,1) double
    options.Log (1,1) logical = false
    options.Step (1,1) double = NaN
end
effectiveHigh=high;
if isfinite(options.Step)
    effectiveHigh=low+floor((high-low)/options.Step+1e-12)*options.Step;
    if effectiveHigh~=high
        warning("radia:optuna:DistributionAdjusted", ...
            "FloatDistribution high was adjusted from %g to %g to align with Step=%g.", ...
            high,effectiveHigh,options.Step);
    end
end
distribution=radia.optuna.internal.DistributionCodec.float( ...
    low,effectiveHigh,options.Log,options.Step);
distribution.name="FloatDistribution";
end
