function study = importStudy(payload, options)
%IMPORTSTUDY camelCase alias of radia.optuna.import_study.
arguments
    payload
    options.Name (1,1) string = ""
    options.StoragePath (1,1) string = ""
    options.Sampler = []
    options.AutoSave (1,1) logical = true
end
study = radia.optuna.import_study(payload, Name=options.Name, ...
    StoragePath=options.StoragePath, Sampler=options.Sampler, ...
    AutoSave=options.AutoSave);
end
