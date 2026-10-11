classdef Boto3ArtifactStore < handle
    %BOTO3ARTIFACTSTORE S3 artifact adapter over an explicit MATLAB client.
    %   The client is a struct of function handles named get_object,
    %   upload_fileobj, and delete_object. No Python fallback is used.

    properties (SetAccess=private)
        bucket (1,1) string
        client
        avoid_buf_copy (1,1) logical
    end

    methods
        function obj=Boto3ArtifactStore(bucket_name,client,options)
            arguments
                bucket_name (1,1) string
                client = []
                options.AvoidBufCopy (1,1) logical = false
            end
            if isempty(client)
                error("radia:optuna:Boto3Dependency", ...
                    "Boto3ArtifactStore requires an explicit MATLAB S3 client adapter.");
            end
            obj.bucket=bucket_name;
            obj.client=client;
            obj.avoid_buf_copy=options.AvoidBufCopy;
        end

        function bytes=open_reader(obj,artifact_id)
            obj.requireOperation("get_object");
            bytes=obj.client.get_object(obj.bucket,string(artifact_id));
            bytes=reshape(uint8(bytes),[],1);
        end

        function write(obj,artifact_id,content_body)
            obj.requireOperation("upload_fileobj");
            bytes=radia.optuna.internal.ArtifactIO.contentBytes(content_body);
            obj.client.upload_fileobj(bytes,obj.bucket,string(artifact_id));
        end

        function remove(obj,artifact_id)
            obj.requireOperation("delete_object");
            obj.client.delete_object(obj.bucket,string(artifact_id));
        end
    end

    methods (Access=private)
        function requireOperation(obj,name)
            if ~isstruct(obj.client) || ~isfield(obj.client,name) || ...
                    ~isa(obj.client.(name),"function_handle")
                error("radia:optuna:Boto3Client", ...
                    "The S3 client adapter must provide %s.",name);
            end
        end
    end
end
