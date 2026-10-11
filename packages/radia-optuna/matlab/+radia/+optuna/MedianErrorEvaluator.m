classdef MedianErrorEvaluator < radia.optuna.BaseErrorEvaluator
    %MEDIANERROREVALUATOR Cache a fraction of the initial median improvement.

    properties (SetAccess=private)
        PairedImprovementEvaluator
        WarmUpTrials (1,1) double
        NInitialTrials (1,1) double
        ThresholdRatio (1,1) double
    end

    properties (Access=private)
        Threshold (1,1) double = NaN
    end

    methods
        function obj=MedianErrorEvaluator(pairedImprovementEvaluator,options)
            arguments
                pairedImprovementEvaluator
                options.WarmUpTrials (1,1) double {mustBeInteger,mustBeNonnegative} = 10
                options.NInitialTrials (1,1) double {mustBeInteger,mustBePositive} = 20
                options.ThresholdRatio (1,1) double {mustBePositive,mustBeFinite} = 0.01
            end
            obj.PairedImprovementEvaluator=pairedImprovementEvaluator;
            obj.WarmUpTrials=options.WarmUpTrials;
            obj.NInitialTrials=options.NInitialTrials;
            obj.ThresholdRatio=options.ThresholdRatio;
        end

        function value=evaluate(obj,trials,study_direction)
            if ~isnan(obj.Threshold)
                value=obj.Threshold;
                return
            end
            if isa(trials,"radia.optuna.Study")
                trials=trials.get_trials();
            end
            if istable(trials)
                trials=trials(trials.State=="COMPLETE",:);
                trials=sortrows(trials,"TrialNumber");
                trialCount=height(trials);
            else
                trials=trials(string({trials.State})=="COMPLETE");
                [~,order]=sort([trials.Number]);
                trials=trials(order);
                trialCount=numel(trials);
            end
            required=obj.WarmUpTrials+obj.NInitialTrials;
            if trialCount<required
                value=-realmin("double");
                return
            end
            criteria=zeros(1,obj.NInitialTrials);
            for index=1:obj.NInitialTrials
                rows=obj.WarmUpTrials+(1:index);
                if istable(trials)
                    selected=trials(rows,:);
                else
                    selected=trials(rows);
                end
                criteria(index)=obj.PairedImprovementEvaluator.evaluate( ...
                    selected,study_direction);
            end
            criteria=sort(criteria);
            medianIndex=floor(numel(criteria)/2)+1;
            obj.Threshold=min(realmax("double"), ...
                criteria(medianIndex)*obj.ThresholdRatio);
            value=obj.Threshold;
        end
    end
end
