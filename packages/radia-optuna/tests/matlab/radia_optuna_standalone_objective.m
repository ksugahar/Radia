function value = radia_optuna_standalone_objective(trial)
%RADIA_OPTUNA_STANDALONE_OBJECTIVE Small wheel-only Simulink block objective.
x = trial.suggestFloat("x", -2, 2);
value = (x - 0.25)^2;
end
