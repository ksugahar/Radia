function study=loadStudy(options)
%LOADSTUDY CamelCase alias for LOAD_STUDY.
arguments
    options.study_name (1,1) string = ""
    options.storage (1,1) string
    options.sampler = []
    options.pruner = []
end
study=radia.optuna.load_study(study_name=options.study_name, ...
    storage=options.storage,sampler=options.sampler,pruner=options.pruner);
end
