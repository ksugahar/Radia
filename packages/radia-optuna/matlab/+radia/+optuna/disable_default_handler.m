function disable_default_handler()
%DISABLE_DEFAULT_HANDLER Disable Optuna's process-local default handler.
radia.optuna.internal.LoggingState.set("default_handler",false);
end
