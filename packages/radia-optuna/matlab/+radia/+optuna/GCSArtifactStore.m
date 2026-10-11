classdef GCSArtifactStore < handle
    %GCSARTIFACTSTORE Google Cloud artifact adapter over a MATLAB client.

    properties (SetAccess=private)
        bucket_name (1,1) string
        client
    end

    methods
        function obj=GCSArtifactStore(bucket_name,client)
            arguments
                bucket_name (1,1) string
                client = []
            end
            if isempty(client)
                error("radia:optuna:GCSDependency", ...
                    "GCSArtifactStore requires an explicit MATLAB GCS client adapter.");
            end
            obj.bucket_name=bucket_name;
            obj.client=client;
        end

        function bytes=open_reader(obj,artifact_id)
            obj.requireOperation("get_blob");
            bytes=obj.client.get_blob(obj.bucket_name,string(artifact_id));
            if isempty(bytes)
                error("radia:optuna:ArtifactNotFound","not found");
            end
            bytes=reshape(uint8(bytes),[],1);
        end

        function write(obj,artifact_id,content_body)
            obj.requireOperation("upload_blob");
            bytes=radia.optuna.internal.ArtifactIO.contentBytes(content_body);
            obj.client.upload_blob(obj.bucket_name,string(artifact_id),bytes);
        end

        function remove(obj,artifact_id)
            obj.requireOperation("delete_blob");
            obj.client.delete_blob(obj.bucket_name,string(artifact_id));
        end
    end

    methods (Access=private)
        function requireOperation(obj,name)
            if ~isstruct(obj.client) || ~isfield(obj.client,name) || ...
                    ~isa(obj.client.(name),"function_handle")
                error("radia:optuna:GCSClient", ...
                    "The GCS client adapter must provide %s.",name);
            end
        end
    end
end
