classdef TrialState < uint8
    %TRIALSTATE Optuna 4.9 trial-state enumeration.

    enumeration
        RUNNING (0)
        COMPLETE (1)
        PRUNED (2)
        FAIL (3)
        WAITING (4)
    end

    methods
        function result=is_finished(obj)
            result=ismember(obj,[radia.optuna.TrialState.COMPLETE, ...
                radia.optuna.TrialState.PRUNED,radia.optuna.TrialState.FAIL]);
        end

        function text=string(obj)
            labels=["RUNNING","COMPLETE","PRUNED","FAIL","WAITING"];
            text=reshape(labels(double(obj)+1),size(obj));
        end

        function text=char(obj)
            if ~isscalar(obj)
                error("radia:optuna:TrialState", ...
                    "char requires a scalar TrialState.");
            end
            text=char(string(obj));
        end
    end

    methods (Static)
        function value=from(value)
            if isa(value,"radia.optuna.TrialState")
                return
            end
            text=upper(string(value));
            valid=["RUNNING","COMPLETE","PRUNED","FAIL","WAITING"];
            if any(~ismember(text,valid))
                error("radia:optuna:TrialState","Unknown trial state.");
            end
            value=repmat(radia.optuna.TrialState.RUNNING,size(text));
            value(text=="COMPLETE")=radia.optuna.TrialState.COMPLETE;
            value(text=="PRUNED")=radia.optuna.TrialState.PRUNED;
            value(text=="FAIL")=radia.optuna.TrialState.FAIL;
            value(text=="WAITING")=radia.optuna.TrialState.WAITING;
        end

        function text=toStorage(value)
            text=string(radia.optuna.TrialState.from(value));
        end
    end
end
