function values = radia_optuna_biobjective(trial)
%RADIA_OPTUNA_BIOBJECTIVE Deterministic two-objective Pareto test surface.
x = trial.suggestFloat("x", -1, 1);
values = [(x + 1)^2, (x - 1)^2];
end
