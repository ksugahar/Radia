classdef FileSystemArtifactStore < handle
    %FILESYSTEMARTIFACTSTORE Store artifact bytes below a checked base path.

    properties (SetAccess=private)
        BasePath (1,1) string
    end

    methods
        function obj=FileSystemArtifactStore(base_path)
            arguments
                base_path (1,1) string
            end
            obj.BasePath=string(base_path);
        end

        function bytes=open_reader(obj,artifact_id)
            path=obj.artifactPath(artifact_id);
            if ~isfile(path)
                error("radia:optuna:ArtifactNotFound","not found");
            end
            bytes=radia.optuna.internal.ArtifactIO.readFile(path);
        end

        function write(obj,artifact_id,content_body)
            path=obj.artifactPath(artifact_id);
            bytes=radia.optuna.internal.ArtifactIO.contentBytes(content_body);
            radia.optuna.internal.ArtifactIO.writeFile(path,bytes,false);
        end

        function remove(obj,artifact_id)
            path=obj.artifactPath(artifact_id);
            if ~isfile(path)
                error("radia:optuna:ArtifactNotFound","not found");
            end
            delete(path);
        end
    end

    methods (Access=private)
        function path=artifactPath(obj,artifact_id)
            artifact_id=string(artifact_id);
            if ~isscalar(artifact_id) || strlength(artifact_id)==0 || ...
                    isfile(artifact_id) || startsWith(artifact_id,["/","\"]) || ...
                    ~isempty(regexp(artifact_id,"^[A-Za-z]:","once")) %#ok<RGXP1>
                error("radia:optuna:ArtifactId", ...
                    "Invalid artifact_id: %s",artifact_id);
            end
            base=string(java.io.File(char(obj.BasePath)).getCanonicalPath());
            candidate=string(java.io.File(char(fullfile(obj.BasePath, ...
                artifact_id))).getCanonicalPath());
            prefix=base+string(filesep);
            if candidate~=base && ~startsWith(candidate,prefix,IgnoreCase=ispc)
                error("radia:optuna:ArtifactId", ...
                    "Invalid artifact_id: %s",artifact_id);
            end
            path=candidate;
        end
    end
end
