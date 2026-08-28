classdef BestValueStagnationEvaluator < radia.optuna.BaseImprovementEvaluator
    %BESTVALUESTAGNATIONEVALUATOR Remaining trials before stagnation limit.

    properties (SetAccess=private)
        MaxStagnationTrials (1,1) double
    end

    methods
        function obj=BestValueStagnationEvaluator(maxStagnationTrials)
            arguments
                maxStagnationTrials (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 30
            end
            obj.MaxStagnationTrials=maxStagnationTrials;
        end

        function remaining=evaluate(obj,trials,direction)
            direction=lower(string(direction));
            if ~ismember(direction,["minimize","maximize"])
                error("radia:optuna:Direction", ...
                    "direction must be minimize or maximize.");
            end
            if istable(trials)
                complete=trials(trials.State=="COMPLETE",:);
                values=complete.Value;
            else
                states=string({trials.State});
                complete=trials(states=="COMPLETE");
                values=[complete.Value]';
            end
            if isempty(values)
                error("radia:optuna:TerminatorTrials", ...
                    "At least one COMPLETE trial is required.");
            end
            if direction=="maximize"
                [~,bestStep]=max(values);
            else
                [~,bestStep]=min(values);
            end
            currentStep=numel(values)-1;
            bestStep=bestStep-1;
            remaining=obj.MaxStagnationTrials-(currentStep-bestStep);
        end
    end
end
