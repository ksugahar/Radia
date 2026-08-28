function summaries=get_all_study_summaries(storage,options)
%GET_ALL_STUDY_SUMMARIES Summarize the zero-or-one study in MAT storage.
arguments
    storage (1,1) string
    options.include_best_trial (1,1) logical = true
end
names=radia.optuna.get_all_study_names(storage);
if isempty(names)
    summaries=radia.optuna.StudySummary.empty(0,1);
    return
end
study=radia.optuna.load_study(study_name=names(1),storage=storage);
best=[];
if options.include_best_trial && isscalar(study.Directions) && ...
        any(study.TrialTable.State=="COMPLETE")
    best=study.best_trial();
end
starts=study.TrialTable.StartTime;
starts=starts(~isnat(starts));
if isempty(starts), started=NaT; else, started=min(starts); end
directions=radia.optuna.StudyDirection.from(study.Directions);
if isscalar(directions)
    summaries=radia.optuna.StudySummary(study.Name,directions,best, ...
        study.UserAttrs,study.SystemAttrs,height(study.TrialTable),started,0);
else
    summaries=radia.optuna.StudySummary(study.Name,[],best, ...
        study.UserAttrs,study.SystemAttrs,height(study.TrialTable),started,0, ...
        directions=directions);
end
end
