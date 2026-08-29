function value = suggestFromDistribution(trial, name, distribution)
%SUGGESTFROMDISTRIBUTION Suggest one value from a typed Optuna distribution.
arguments
    trial (1,1) radia.optuna.Trial
    name (1,1) string
    distribution
end

switch distribution.kind
    case "categorical"
        value = trial.suggest_categorical(name, distribution.choices);
    case "integer"
        value = trial.suggest_int(name, distribution.low, ...
            distribution.high, log=distribution.log, ...
            step=distribution.step);
    otherwise
        value = trial.suggest_float(name, distribution.low, ...
            distribution.high, log=distribution.log, ...
            step=distribution.step);
end
end
