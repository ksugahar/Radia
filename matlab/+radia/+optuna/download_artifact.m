function download_artifact(options)
%DOWNLOAD_ARTIFACT Copy an artifact to a new local file.
arguments
    options.artifact_store
    options.file_path (1,1) string
    options.artifact_id (1,1) string
end
bytes=options.artifact_store.open_reader(options.artifact_id);
radia.optuna.internal.ArtifactIO.writeFile( ...
    options.file_path,bytes,true);
end
