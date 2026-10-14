function disable_propagation()
%DISABLE_PROPAGATION Disable propagation from the root Optuna logger.
radia.optuna.internal.LoggingState.set("propagation",false);
end
