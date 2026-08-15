function distribution=CategoricalDistribution(choices)
%CATEGORICALDISTRIBUTION Public Optuna-style categorical descriptor.
distribution=radia.optuna.internal.DistributionCodec.categorical(choices);
end
