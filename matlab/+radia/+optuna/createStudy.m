function study = createStudy(options)
%CREATESTUDY CamelCase alias for CREATE_STUDY.
arguments
    options.study_name (1,1) string = "radia-study"
    options.direction (1,1) string = "minimize"
    options.directions string = string.empty
    options.sampler = []
    options.pruner = []
    options.storage (1,1) string = ""
    options.StoragePath (1,1) string = ""
    options.AutoSave (1,1) logical = true
    options.ProgressFcn = []
end
study = radia.optuna.create_study( ...
    study_name=options.study_name, direction=options.direction, ...
    directions=options.directions, sampler=options.sampler, pruner=options.pruner, ...
    storage=options.storage, StoragePath=options.StoragePath, AutoSave=options.AutoSave, ...
    ProgressFcn=options.ProgressFcn);
end
