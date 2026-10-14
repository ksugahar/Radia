function deleteStudy(options)
%DELETESTUDY CamelCase alias for DELETE_STUDY.
arguments
    options.study_name (1,1) string
    options.storage (1,1) string
end
radia.optuna.delete_study(study_name=options.study_name,storage=options.storage);
end
