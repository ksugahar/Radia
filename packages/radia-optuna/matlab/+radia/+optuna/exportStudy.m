function payload = exportStudy(study, options)
%EXPORTSTUDY camelCase alias of radia.optuna.export_study.
arguments
    study (1,1) radia.optuna.Study
    options.Path (1,1) string = ""
end
payload = radia.optuna.export_study(study, Path=options.Path);
end
