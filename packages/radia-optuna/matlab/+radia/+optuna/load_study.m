function study=load_study(options)
%LOAD_STUDY Load one MAT-backed Optuna-compatible study.
arguments
    options.study_name (1,1) string = ""
    options.storage (1,1) string
    options.sampler = []
    options.pruner = []
end
if ~isfile(options.storage) && ~isfile(options.storage+".bak")
    error("radia:optuna:StudyNotFound", ...
        "Study storage '%s' does not exist.",options.storage);
end
study=radia.optuna.Study(Name=options.study_name, ...
    StoragePath=options.storage,Sampler=options.sampler, ...
    Pruner=options.pruner);
if strlength(options.study_name)>0 && study.Name~=options.study_name
    error("radia:optuna:StudyNotFound", ...
        "Storage contains study '%s', not '%s'.",study.Name,options.study_name);
end
end
