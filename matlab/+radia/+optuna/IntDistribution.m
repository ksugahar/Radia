function distribution=IntDistribution(low,high,options)
%INTDISTRIBUTION Public Optuna-style integer distribution descriptor.
arguments
    low (1,1) double
    high (1,1) double
    options.Log (1,1) logical = false
    options.Step (1,1) double = 1
end
if low~=floor(low) || high~=floor(high) || ...
        options.Step~=floor(options.Step) || options.Step<=0
    error("radia:optuna:IntegerDistribution", ...
        "Integer bounds and Step must be integers and Step must be positive.");
end
if options.Log && options.Step~=1
    error("radia:optuna:IntegerDistribution", ...
        "Log integer distributions require Step=1.");
end
effectiveHigh=low+floor((high-low)/options.Step)*options.Step;
distribution=radia.optuna.internal.DistributionCodec.integer( ...
    low,effectiveHigh,options.Log,options.Step);
end
