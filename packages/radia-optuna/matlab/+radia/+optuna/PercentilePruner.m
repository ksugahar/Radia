classdef PercentilePruner < radia.optuna.BasePruner
    %PERCENTILEPRUNER Keep the requested percentile of completed trials.

    properties (SetAccess=private)
        Percentile (1,1) double
        NStartupTrials (1,1) double
        NWarmupSteps (1,1) double
        IntervalSteps (1,1) double
        NMinTrials (1,1) double
    end

    methods
        function obj=PercentilePruner(percentile,options)
            arguments
                percentile (1,1) double
                options.NStartupTrials (1,1) double = 5
                options.NWarmupSteps (1,1) double = 0
                options.IntervalSteps (1,1) double = 1
                options.NMinTrials (1,1) double = 1
            end
            if percentile<0 || percentile>100 || ~isfinite(percentile)
                error("radia:optuna:PrunerPercentile", ...
                    "Percentile must be between zero and 100 inclusive.");
            end
            if options.NStartupTrials<0 || ...
                    options.NStartupTrials~=floor(options.NStartupTrials)
                error("radia:optuna:PrunerStartup", ...
                    "NStartupTrials must be a nonnegative integer.");
            end
            radia.optuna.internal.PrunerSupport.validateSchedule( ...
                options.NWarmupSteps,options.IntervalSteps);
            if options.NMinTrials<1 || ...
                    options.NMinTrials~=floor(options.NMinTrials)
                error("radia:optuna:PrunerMinimum", ...
                    "NMinTrials must be a positive integer.");
            end
            obj.Percentile=percentile;
            obj.NStartupTrials=options.NStartupTrials;
            obj.NWarmupSteps=options.NWarmupSteps;
            obj.IntervalSteps=options.IntervalSteps;
            obj.NMinTrials=options.NMinTrials;
        end

        function decision=shouldPrune(obj,study,trial)
            decision=false;
            completed=study.TrialTable.State=="COMPLETE";
            if sum(completed)<obj.NStartupTrials || ~any(completed) || ...
                    isempty(trial.IntermediateValues)
                return
            end
            if ~radia.optuna.internal.PrunerSupport.isFirstInInterval( ...
                    trial,obj.NWarmupSteps,obj.IntervalSteps)
                return
            end
            best=radia.optuna.internal.PrunerSupport.bestIntermediate( ...
                trial,study.Directions(1));
            if isnan(best)
                decision=true;
                return
            end
            step=max(trial.IntermediateValues.Step);
            reference=study.intermediateValuesAtStep(step);
            if numel(reference)<obj.NMinTrials
                return
            end
            percentile=obj.Percentile;
            if study.Directions(1)=="maximize"
                percentile=100-percentile;
            end
            threshold=radia.optuna.internal.PrunerSupport.percentile( ...
                reference,percentile);
            if isnan(threshold)
                return
            end
            if study.Directions(1)=="maximize"
                decision=best<threshold;
            else
                decision=best>threshold;
            end
        end
    end
end
