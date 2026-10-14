function value = radia_optuna_session_quadratic(values)
%RADIA_OPTUNA_SESSION_QUADRATIC Model-shaped objective for block tests.
value = (values.simulink_session_x - 0.2)^2;
end
