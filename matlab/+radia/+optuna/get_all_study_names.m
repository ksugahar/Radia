function names=get_all_study_names(storage)
%GET_ALL_STUDY_NAMES Return the zero-or-one study name in a MAT storage.
arguments
    storage (1,1) string
end
if ~isfile(storage) && ~isfile(storage+".bak")
    names=strings(0,1);
    return
end
study=radia.optuna.Study(StoragePath=storage,AutoSave=false);
names=study.Name;
end
