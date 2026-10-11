function copyStudy(options)
%COPYSTUDY CamelCase alias for COPY_STUDY.
arguments
    options.from_study_name (1,1) string
    options.from_storage (1,1) string
    options.to_storage (1,1) string
    options.to_study_name (1,1) string = ""
end
radia.optuna.copy_study(from_study_name=options.from_study_name, ...
    from_storage=options.from_storage,to_storage=options.to_storage, ...
    to_study_name=options.to_study_name);
end
