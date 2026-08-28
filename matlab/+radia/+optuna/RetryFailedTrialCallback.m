classdef RetryFailedTrialCallback < radia.optuna.RetryHeartbeatStaleTrialCallback
    %RETRYFAILEDTRIALCALLBACK Deprecated Optuna 4.9 retry alias.

    methods
        function obj=RetryFailedTrialCallback(options)
            arguments
                options.max_retry double = NaN
                options.inherit_intermediate_values (1,1) logical = false
            end
            warningState=warning("off","radia:optuna:ExperimentalWarning");
            cleanup=onCleanup(@()warning(warningState));
            obj@radia.optuna.RetryHeartbeatStaleTrialCallback( ...
                max_retry=options.max_retry, ...
                inherit_intermediate_values= ...
                options.inherit_intermediate_values);
            clear cleanup
            warning("radia:optuna:FutureWarning", ...
                "RetryFailedTrialCallback is deprecated; use " + ...
                "RetryHeartbeatStaleTrialCallback instead.");
        end
    end
end
