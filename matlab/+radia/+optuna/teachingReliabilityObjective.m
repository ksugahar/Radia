function result = teachingReliabilityObjective(values, context)
%TEACHINGRELIABILITYOBJECTIVE Deterministic complete/pruned/failed exercise.
%   A brute-force grid over [-1,-0.5,0,0.5,1] produces at least one failed
%   trial, one pruned trial, and successful trials without relying on timing.
arguments
    values (1,1) struct
    context (1,1) radia.optuna.ObjectiveContext
end
x = double(values.x);
metric = (x - 0.1)^2;
context.report(metric + 0.1, 0);
if x < -0.5
    error("radia:optuna:TeachingFailure", ...
        "Deliberate teaching failure for x=%g.", x);
end
if x > 0.5
    context.prune();
end
context.report(metric, 1);
result = struct( ...
    "Objectives", metric, ...
    "InequalityConstraints", x - 0.5, ...
    "UserData", struct("exercise", "prune-and-fail"));
end
