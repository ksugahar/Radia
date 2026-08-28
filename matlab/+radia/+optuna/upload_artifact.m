function artifact_id=upload_artifact(options)
%UPLOAD_ARTIFACT Store a file and attach its metadata to a study or trial.
arguments
    options.artifact_store
    options.file_path (1,1) string
    options.study_or_trial
    options.storage = []
    options.mimetype (1,1) string = ""
    options.encoding (1,1) string = ""
end
if ~isfile(options.file_path)
    error("radia:optuna:ArtifactRead", ...
        "Artifact file does not exist: %s",options.file_path);
end
[~,name,extension]=fileparts(options.file_path);
filename=name+extension;
artifact_id=string(java.util.UUID.randomUUID());
mimetype=options.mimetype;
encoding=options.encoding;
if strlength(mimetype)==0
    mimetype=guessMimeType(extension);
end
if strlength(encoding)==0
    encoding=guessEncoding(extension);
end
meta=radia.optuna.ArtifactMeta(artifact_id,filename,mimetype,encoding);
payload=struct("artifact_id",meta.artifact_id, ...
    "filename",meta.filename,"mimetype",meta.mimetype, ...
    "encoding",meta.encoding);
key="artifacts:"+artifact_id;
target=options.study_or_trial;
if isa(target,"radia.optuna.FrozenTrial")
    if isempty(options.storage)
        error("radia:optuna:ArtifactStorage", ...
            "storage is required for FrozenTrial.");
    end
    options.storage.set_trial_system_attr(target.Number,key,jsonencode(payload));
elseif ismethod(target,"set_system_attr")
    target.set_system_attr(key,jsonencode(payload));
else
    error("radia:optuna:ArtifactTarget", ...
        "study_or_trial must be an Optuna Study, Trial, or FrozenTrial.");
end
bytes=radia.optuna.internal.ArtifactIO.readFile(options.file_path);
options.artifact_store.write(artifact_id,bytes);
end

function value=guessMimeType(extension)
switch lower(extension)
    case ".txt", value="text/plain";
    case ".json", value="application/json";
    case ".csv", value="text/csv";
    case ".pdf", value="application/pdf";
    case ".png", value="image/png";
    case {".jpg",".jpeg"}, value="image/jpeg";
    case ".gz", value="application/octet-stream";
    otherwise, value="application/octet-stream";
end
end

function value=guessEncoding(extension)
if lower(extension)==".gz"
    value="gzip";
else
    value=missing;
end
end
