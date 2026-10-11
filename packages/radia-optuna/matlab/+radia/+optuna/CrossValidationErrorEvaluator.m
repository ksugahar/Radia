classdef CrossValidationErrorEvaluator < radia.optuna.BaseErrorEvaluator
    %CROSSVALIDATIONERROREVALUATOR Scaled variance of the best trial's folds.

    methods
        function value=evaluate(~,trials,study_direction)
            complete=radia.optuna.CrossValidationErrorEvaluator. ...
                completeTrials(trials);
            if isempty(complete)
                error("radia:optuna:TerminatorTrials", ...
                    "At least one COMPLETE trial is required.");
            end
            values=reshape([complete.Value],[],1);
            direction=radia.optuna.StudyDirection.from(study_direction);
            if direction==radia.optuna.StudyDirection.MAXIMIZE
                [~,index]=max(values);
            else
                [~,index]=min(values);
            end
            key=matlab.lang.makeValidName("terminator:cv_scores");
            attributes=complete(index).SystemAttrs;
            if ~isfield(attributes,key)
                error("radia:optuna:CrossValidationScores", ...
                    "Cross-validation scores have not been reported. Call report_cross_validation_scores during the trial.");
            end
            scores=reshape(double(attributes.(key)),1,[]);
            if numel(scores)<=1
                error("radia:optuna:CrossValidationScores", ...
                    "The length of scores must be greater than one.");
            end
            count=numel(scores);
            value=sqrt((1/count+1/(count-1))*var(scores,1));
        end
    end

    methods (Static, Access=private)
        function complete=completeTrials(trials)
            if isa(trials,"radia.optuna.Study")
                trials=trials.get_trials();
            end
            if istable(trials)
                error("radia:optuna:TerminatorTrials", ...
                    "Cross-validation evaluation requires FrozenTrial values with system attributes.");
            end
            if isempty(trials)
                complete=trials;
                return
            end
            complete=trials(string({trials.State})=="COMPLETE");
        end
    end
end
