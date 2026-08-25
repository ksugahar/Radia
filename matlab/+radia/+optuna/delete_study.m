function delete_study(options)
%DELETE_STUDY Delete one MAT-backed study and its recovery backup.
arguments
    options.study_name (1,1) string
    options.storage (1,1) string
end
study=radia.optuna.load_study(study_name=options.study_name, ...
    storage=options.storage);
path=study.StoragePath;
if isfile(path), delete(path); end
if isfile(path+".bak"), delete(path+".bak"); end
end
