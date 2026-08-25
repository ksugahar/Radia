function study = create_study(options)
%CREATE_STUDY Create a study from the verified Optuna 4.9 MATLAB subset.
arguments
    options.study_name (1,1) string = ""
    options.direction = string.empty
    options.directions = string.empty
    options.sampler = []
    options.pruner = []
    options.storage (1,1) string = ""
    options.StoragePath (1,1) string = ""
    options.AutoSave (1,1) logical = true
    options.ProgressFcn = []
    options.load_if_exists (1,1) logical = false
end
if ~isempty(options.directions) && ~isempty(options.direction)
    error("radia:optuna:Direction", ...
        "Specify only one of direction or directions.");
elseif ~isempty(options.directions)
    directions = options.directions;
elseif ~isempty(options.direction)
    directions = options.direction;
else
    directions = "minimize";
end
storage = options.storage;
if strlength(options.StoragePath) > 0
    storage = options.StoragePath;
end
if strlength(storage)>0 && (isfile(storage) || isfile(storage+".bak"))
    if ~options.load_if_exists
        error("radia:optuna:DuplicatedStudy", ...
            "Study storage '%s' already exists. Set load_if_exists=true to load it.", ...
            storage);
    end
    study=radia.optuna.load_study(study_name=options.study_name, ...
        storage=storage,sampler=options.sampler,pruner=options.pruner);
    return
end
study = radia.optuna.Study( ...
    Name=options.study_name, Directions=directions, Sampler=options.sampler, ...
    Pruner=options.pruner, StoragePath=storage, AutoSave=options.AutoSave, ...
    ProgressFcn=options.ProgressFcn);
end
