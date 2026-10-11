function result = teachingParetoObjective(values, context)
%TEACHINGPARETOOBJECTIVE Two competing objectives with a visible Pareto set.
%   For 0 <= x <= 1, the objectives x^2 and (x-1)^2 trade off directly.
arguments
    values (1,1) struct
    context (1,1) radia.optuna.ObjectiveContext
end
x = double(values.x);
result = struct( ...
    "Objectives", [x^2, (x - 1)^2], ...
    "InequalityConstraints", zeros(1,0), ...
    "UserData", struct("exercise", "pareto-selection"));
context.setUserAttr("teaching_exercise", "pareto-selection");
end
