function result=get_all_artifact_meta(study_or_trial,options)
%GET_ALL_ARTIFACT_META Return metadata attached to a study or trial.
arguments
    study_or_trial
    options.storage = []
end
if isa(study_or_trial,"radia.optuna.FrozenTrial")
    if isempty(options.storage)
        error("radia:optuna:ArtifactStorage", ...
            "storage is required for FrozenTrial.");
    end
    attributes=options.storage.get_trial_system_attrs( ...
        study_or_trial.Number);
elseif ismethod(study_or_trial,"system_attrs")
    attributes=study_or_trial.system_attrs();
else
    error("radia:optuna:ArtifactTarget", ...
        "study_or_trial must be an Optuna Study, Trial, or FrozenTrial.");
end
names=string(fieldnames(attributes));
result=radia.optuna.ArtifactMeta.empty(0,1);
for index=1:numel(names)
    original=names(index);
    if ~startsWith(original,"artifacts_") && ...
            ~startsWith(original,"artifacts:")
        continue
    end
    payload=jsondecode(char(string(attributes.(original))));
    encoding=missing;
    if isfield(payload,"encoding") && ~isempty(payload.encoding)
        encoding=string(payload.encoding);
    end
    result(end+1,1)=radia.optuna.ArtifactMeta( ...
        string(payload.artifact_id),string(payload.filename), ...
        string(payload.mimetype),encoding); %#ok<AGROW>
end
end
