function distribution=FloatDistribution(low,high,options)
%FLOATDISTRIBUTION Public Optuna-style float distribution descriptor.
arguments
    low (1,1) double
    high (1,1) double
    options.Log (1,1) logical = false
    options.Step (1,1) double = NaN
end
distribution=radia.optuna.internal.DistributionCodec.float( ...
    low,high,options.Log,options.Step);
end
