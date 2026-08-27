function set_verbosity(verbosity)
%SET_VERBOSITY Set the process-local Optuna logging level.
arguments
    verbosity (1,1) double
end
radia.optuna.internal.LoggingState.set("verbosity",verbosity);
end
