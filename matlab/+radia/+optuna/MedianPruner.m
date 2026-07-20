classdef MedianPruner < handle
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
            latestStep = trial.IntermediateValues.Step(end);
            latestValue = trial.IntermediateValues.Value(end);
            if latestStep < obj.NWarmupSteps || ...
                    mod(latestStep - obj.NWarmupSteps, obj.IntervalSteps) ~= 0
                return;
            end
            completed = study.TrialTable.State == "COMPLETE";
            if sum(completed) < max(obj.NStartupTrials, obj.MinCompletedTrials)
                return;
            end
            reference = study.intermediateValuesAtStep(latestStep);
            if numel(reference) < obj.MinCompletedTrials || isempty(reference)
                return;
            end
            threshold = median(reference, "omitnan");
            if study.Directions(1) == "minimize"
                decision = latestValue > threshold;
            else
                decision = latestValue < threshold;
            end
        end
    end
end
