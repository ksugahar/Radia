function distribution=json_to_distribution(encoded)
%JSON_TO_DISTRIBUTION Deserialize Optuna 4.9 distribution JSON.
arguments
    encoded {mustBeTextScalar}
end
try
    payload=jsondecode(char(encoded));
catch cause
    error("radia:optuna:DistributionJSON", ...
        "Invalid distribution JSON: %s",cause.message);
end
if ~isstruct(payload) || ~all(isfield(payload,["name","attributes"]))
    error("radia:optuna:DistributionJSON", ...
        "Distribution JSON requires name and attributes.");
end
name=string(payload.name);
attributes=payload.attributes;
switch name
    case "FloatDistribution"
        step=NaN;
        if ~isempty(attributes.step), step=double(attributes.step); end
        distribution=radia.optuna.FloatDistribution( ...
            attributes.low,attributes.high,Log=attributes.log,Step=step);
    case "IntDistribution"
        distribution=radia.optuna.IntDistribution( ...
            attributes.low,attributes.high,Log=attributes.log, ...
            Step=attributes.step);
    case "CategoricalDistribution"
        distribution=radia.optuna.CategoricalDistribution(attributes.choices);
    case "UniformDistribution"
        distribution=radia.optuna.FloatDistribution( ...
            attributes.low,attributes.high);
    case "LogUniformDistribution"
        distribution=radia.optuna.FloatDistribution( ...
            attributes.low,attributes.high,Log=true);
    case "DiscreteUniformDistribution"
        distribution=radia.optuna.FloatDistribution( ...
            attributes.low,attributes.high,Step=attributes.q);
    case "IntUniformDistribution"
        distribution=radia.optuna.IntDistribution( ...
            attributes.low,attributes.high,Step=attributes.step);
    case "IntLogUniformDistribution"
        distribution=radia.optuna.IntDistribution( ...
            attributes.low,attributes.high,Log=true,Step=attributes.step);
    otherwise
        error("radia:optuna:DistributionJSON", ...
            "Unknown distribution '%s'.",name);
end
distribution.name=name;
end
