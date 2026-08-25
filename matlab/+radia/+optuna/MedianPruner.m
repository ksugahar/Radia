classdef MedianPruner < radia.optuna.BasePruner
    %MEDIANPRUNER Stop trials that are worse than completed-trial median.

    properties
        NStartupTrials (1,1) double = 5
        NWarmupSteps (1,1) double = 0
        IntervalSteps (1,1) double = 1
        MinCompletedTrials (1,1) double = 1
    end

    methods
        function obj = MedianPruner(options)
            arguments
                options.NStartupTrials (1,1) double = 5
                options.NWarmupSteps (1,1) double = 0
                options.IntervalSteps (1,1) double = 1
                options.MinCompletedTrials (1,1) double = 1
            end
            if options.NStartupTrials < 0 || ...
                    options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:PrunerStartup", ...
                    "NStartupTrials must be a nonnegative integer.");
            end
            radia.optuna.internal.PrunerSupport.validateSchedule( ...
                options.NWarmupSteps, options.IntervalSteps);
            if options.MinCompletedTrials < 1 || ...
                    options.MinCompletedTrials ~= floor(options.MinCompletedTrials)
                error("radia:optuna:PrunerMinimum", ...
                    "MinCompletedTrials must be a positive integer.");
            end
            obj.NStartupTrials = options.NStartupTrials;
            obj.NWarmupSteps = options.NWarmupSteps;
            obj.IntervalSteps = options.IntervalSteps;
            obj.MinCompletedTrials = options.MinCompletedTrials;
        end

        function decision = shouldPrune(obj, study, trial)
            decision = false;
            if trial.State ~= "RUNNING" || isempty(trial.IntermediateValues)
                return;
            end
            if numel(study.Directions) ~= 1
                error("radia:optuna:PrunerMultiObjective", ...
                    "MedianPruner supports single-objective studies only.");
            end
            if ~radia.optuna.internal.PrunerSupport.isFirstInInterval( ...
                    trial, obj.NWarmupSteps, obj.IntervalSteps)
                return;
            end
            latestStep = max(trial.IntermediateValues.Step);
            completed = study.TrialTable.State == "COMPLETE";
            if sum(completed) < max(obj.NStartupTrials, obj.MinCompletedTrials)
                return;
            end
            reference = study.intermediateValuesAtStep(latestStep);
            if numel(reference) < obj.MinCompletedTrials || isempty(reference)
                return;
            end
            threshold = radia.optuna.internal.PrunerSupport.percentile( ...
                reference, 50);
            bestIntermediate = ...
                radia.optuna.internal.PrunerSupport.bestIntermediate( ...
                trial, study.Directions(1));
            if isnan(bestIntermediate)
                decision = true;
                return
            end
            if isnan(threshold)
                return
            end
            if study.Directions(1) == "minimize"
                decision = bestIntermediate > threshold;
            else
                decision = bestIntermediate < threshold;
            end
        end
    end
end
