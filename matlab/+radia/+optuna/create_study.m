function study = create_study(options)
%CREATE_STUDY Optuna-compatible MATLAB study factory.
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
if ~isempty(options.directions)
    directions = options.directions;
else
    directions = options.direction;
end
storage = options.storage;
if strlength(options.StoragePath) > 0
    storage = options.StoragePath;
end
study = radia.optuna.Study( ...
    Name=options.study_name, Directions=directions, Sampler=options.sampler, ...
    Pruner=options.pruner, StoragePath=storage, AutoSave=options.AutoSave, ...
    ProgressFcn=options.ProgressFcn);
end
