function value=radia_optuna_quadratic(trial)
%RADIA_OPTUNA_QUADRATIC Deterministic test objective for the Simulink block.
x=trial.suggestFloat("x",-4,4); value=(x-1.25)^2;
end
