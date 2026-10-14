classdef ArtifactMeta
    %ARTIFACTMETA Metadata associated with one uploaded artifact.

    properties
        artifact_id (1,1) string
        filename (1,1) string
        mimetype (1,1) string
        encoding (1,1) string
    end

    methods
        function obj=ArtifactMeta(artifact_id,filename,mimetype,encoding)
            arguments
                artifact_id (1,1) string
                filename (1,1) string
                mimetype (1,1) string
                encoding (1,1) string = missing
            end
            obj.artifact_id=artifact_id;
            obj.filename=filename;
            obj.mimetype=mimetype;
            obj.encoding=encoding;
        end
    end
end
