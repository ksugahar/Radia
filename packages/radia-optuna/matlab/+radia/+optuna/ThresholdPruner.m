classdef ThresholdPruner < radia.optuna.BasePruner
    %THRESHOLDPRUNER Prune when the latest value crosses fixed bounds.

    properties (SetAccess=private)
        Lower (1,1) double
        Upper (1,1) double
        NWarmupSteps (1,1) double
        IntervalSteps (1,1) double
    end

    methods
        function obj=ThresholdPruner(options)
            arguments
                options.Lower (1,1) double = NaN
                options.Upper (1,1) double = NaN
                options.NWarmupSteps (1,1) double = 0
                options.IntervalSteps (1,1) double = 1
            end
            if isnan(options.Lower) && isnan(options.Upper)
                error("radia:optuna:PrunerThreshold", ...
                    "Either Lower or Upper must be specified.");
            end
            lower=options.Lower;
            upper=options.Upper;
            if isnan(lower), lower=-Inf; end
            if isnan(upper), upper=Inf; end
            if lower>upper
                error("radia:optuna:PrunerThreshold", ...
                    "Lower and Upper must be ordered non-NaN thresholds.");
            end
            radia.optuna.internal.PrunerSupport.validateSchedule( ...
                options.NWarmupSteps,options.IntervalSteps);
            obj.Lower=lower;
            obj.Upper=upper;
            obj.NWarmupSteps=options.NWarmupSteps;
            obj.IntervalSteps=options.IntervalSteps;
        end

        function decision=shouldPrune(obj,~,trial)
            decision=false;
            if isempty(trial.IntermediateValues) || ...
                    ~radia.optuna.internal.PrunerSupport.isFirstInInterval( ...
                    trial,obj.NWarmupSteps,obj.IntervalSteps)
                return
            end
            value=radia.optuna.internal.PrunerSupport.latestIntermediate(trial);
            decision=isnan(value) || value<obj.Lower || value>obj.Upper;
        end
    end
end
