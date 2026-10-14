function pruner = prunerFromName(name)
%PRUNERFROMNAME Construct a pruner from its public short name.
arguments
    name (1,1) string
end

switch lower(name)
    case "none"
        pruner = [];
    case "median"
        pruner = radia.optuna.MedianPruner();
    case "hyperband"
        pruner = radia.optuna.HyperbandPruner();
    case "percentile"
        pruner = radia.optuna.PercentilePruner(25);
    case "patient"
        pruner = radia.optuna.PatientPruner( ...
            radia.optuna.MedianPruner(), 1);
    case "successivehalving"
        pruner = radia.optuna.SuccessiveHalvingPruner();
    case "threshold"
        pruner = radia.optuna.ThresholdPruner();
    otherwise
        error("radia:optuna:PrunerName", ...
            "Unknown pruner '%s'.", name);
end
end
