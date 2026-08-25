function summaries=getAllStudySummaries(storage,options)
%GETALLSTUDYSUMMARIES CamelCase alias for GET_ALL_STUDY_SUMMARIES.
arguments
    storage (1,1) string
    options.IncludeBestTrial (1,1) logical = true
end
summaries=radia.optuna.get_all_study_summaries(storage, ...
    include_best_trial=options.IncludeBestTrial);
end
