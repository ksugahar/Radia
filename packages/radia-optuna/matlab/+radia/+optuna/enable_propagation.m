function enable_propagation()
%ENABLE_PROPAGATION Enable propagation from the root Optuna logger.
radia.optuna.internal.LoggingState.set("propagation",true);
end
