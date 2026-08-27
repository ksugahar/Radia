function encoded=distribution_to_json(distribution)
%DISTRIBUTION_TO_JSON Serialize a distribution using Optuna 4.9 JSON.
arguments
    distribution
end
if ~radia.optuna.internal.DistributionCodec.isSpec(distribution)
    error("radia:optuna:Distribution", ...
        "A distribution requires kind and name fields.");
end
distribution=radia.optuna.internal.DistributionCodec.normalize(distribution);
name=string(distribution.name);
switch name
    case "FloatDistribution"
        attributes=struct("step",distribution.step,"low",distribution.low, ...
            "high",distribution.high,"log",distribution.log);
    case "IntDistribution"
        attributes=struct("log",distribution.log,"step",distribution.step, ...
            "low",distribution.low,"high",distribution.high);
    case "CategoricalDistribution"
        attributes=struct("choices",{distribution.choices});
    case {"UniformDistribution","LogUniformDistribution"}
        attributes=struct("low",distribution.low,"high",distribution.high);
    case "DiscreteUniformDistribution"
        attributes=struct("low",distribution.low,"high",distribution.high, ...
            "q",distribution.step);
    case {"IntUniformDistribution","IntLogUniformDistribution"}
        attributes=struct("step",distribution.step,"low",distribution.low, ...
            "high",distribution.high);
    otherwise
        error("radia:optuna:DistributionKind", ...
            "Unsupported distribution '%s'.",name);
end
encoded=string(jsonencode(struct("name",name,"attributes",attributes)));
end
