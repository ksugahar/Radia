function searchSpace=intersection_search_space(trials,options)
%INTERSECTION_SEARCH_SPACE Return distributions common to finished trials.
arguments
    trials
    options.include_pruned (1,1) logical = false
end
searchSpace=radia.optuna.IntersectionSearchSpace.fromTrials( ...
    trials,options.include_pruned);
end
