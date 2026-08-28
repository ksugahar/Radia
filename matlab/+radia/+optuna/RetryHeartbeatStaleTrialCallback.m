classdef RetryHeartbeatStaleTrialCallback
    %RETRYHEARTBEATSTALETRIALCALLBACK Requeue heartbeat-stale trials.

    properties (SetAccess=immutable)
        max_retry double
        inherit_intermediate_values (1,1) logical
    end

    methods
        function obj=RetryHeartbeatStaleTrialCallback(options)
            arguments
                options.max_retry double = NaN
                options.inherit_intermediate_values (1,1) logical = false
            end
            if ~isnan(options.max_retry) && ...
                    (~isscalar(options.max_retry) || ...
                    options.max_retry~=floor(options.max_retry))
                error("radia:optuna:RetryCount", ...
                    "max_retry must be an integer or NaN for no limit.");
            end
            warning("radia:optuna:ExperimentalWarning", ...
                "RetryHeartbeatStaleTrialCallback is experimental.");
            obj.max_retry=options.max_retry;
            obj.inherit_intermediate_values= ...
                options.inherit_intermediate_values;
        end

        function call(obj,study,trial)
            arguments
                obj
                study (1,1) radia.optuna.Study
                trial (1,1) radia.optuna.FrozenTrial
            end
            systemAttrs=trial.SystemAttrs;
            if ~isfield(systemAttrs,"failed_trial")
                systemAttrs.failed_trial=trial.Number;
            end
            if ~isfield(systemAttrs,"retry_history")
                systemAttrs.retry_history=zeros(1,0);
                sharedHistory=false;
            else
                sharedHistory=true;
            end
            history=reshape(double(systemAttrs.retry_history),1,[]);
            history(end+1)=trial.Number;
            systemAttrs.retry_history=history;
            if sharedHistory
                % Python's shallow dict expansion retains the history list,
                % so appending is observable on the supplied FrozenTrial.
                trial.setSystemAttr("retry_history",history);
            end
            if ~isnan(obj.max_retry) && obj.max_retry<numel(history)
                return
            end
            intermediate=radia.optuna.Trial.emptyIntermediateTable();
            if obj.inherit_intermediate_values
                intermediate=trial.IntermediateValues;
            end
            study.add_trial(radia.optuna.create_trial( ...
                state="WAITING",params=trial.Params, ...
                distributions=trial.Distributions, ...
                intermediate_values=intermediate, ...
                user_attrs=trial.UserAttrs,system_attrs=systemAttrs));
        end
    end

    methods (Static)
        function number=retried_trial_number(trial)
            warning("radia:optuna:ExperimentalWarning", ...
                "retried_trial_number is experimental.");
            if isfield(trial.SystemAttrs,"failed_trial")
                number=double(trial.SystemAttrs.failed_trial);
            else
                number=[];
            end
        end

        function history=retry_history(trial)
            warning("radia:optuna:ExperimentalWarning", ...
                "retry_history is experimental.");
            if isfield(trial.SystemAttrs,"retry_history")
                history=reshape(double(trial.SystemAttrs.retry_history),1,[]);
            else
                history=zeros(1,0);
            end
        end
    end
end
