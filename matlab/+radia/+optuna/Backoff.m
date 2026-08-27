classdef Backoff < handle
    %BACKOFF Exponential-retry middleware for an artifact store.

    properties (SetAccess=private)
        Backend
        MaxRetries (1,1) double
        Multiplier (1,1) double
        MinDelay (1,1) double
        MaxDelay (1,1) double
    end

    methods
        function obj=Backoff(backend,options)
            arguments
                backend
                options.MaxRetries (1,1) double {mustBeInteger,mustBePositive} = 10
                options.Multiplier (1,1) double {mustBePositive} = 2
                options.MinDelay (1,1) double {mustBePositive} = 0.1
                options.MaxDelay (1,1) double {mustBePositive} = 30
            end
            if options.MaxDelay<=options.MinDelay
                error("radia:optuna:ArtifactBackoff", ...
                    "MaxDelay must be greater than MinDelay.");
            end
            obj.Backend=backend;
            obj.MaxRetries=options.MaxRetries;
            obj.Multiplier=options.Multiplier;
            obj.MinDelay=options.MinDelay;
            obj.MaxDelay=options.MaxDelay;
        end

        function bytes=open_reader(obj,artifact_id)
            bytes=obj.retry("open_reader",artifact_id,[]);
        end

        function write(obj,artifact_id,content_body)
            obj.retry("write",artifact_id,content_body);
        end

        function remove(obj,artifact_id)
            % Optuna 4.9 calls remove on every retry slot even after success.
            for index=0:obj.MaxRetries-1
                try
                    obj.Backend.remove(artifact_id);
                catch exception
                    if string(exception.identifier)=="radia:optuna:ArtifactNotFound" || ...
                            index==obj.MaxRetries-1
                        rethrow(exception);
                    end
                    pause(obj.sleepSeconds(index));
                end
            end
        end
    end

    methods (Access=private)
        function output=retry(obj,operation,artifact_id,content)
            output=[];
            for index=0:obj.MaxRetries-1
                try
                    if operation=="open_reader"
                        output=obj.Backend.open_reader(artifact_id);
                    else
                        obj.Backend.write(artifact_id,content);
                    end
                    return
                catch exception
                    if string(exception.identifier)=="radia:optuna:ArtifactNotFound" || ...
                            index==obj.MaxRetries-1
                        rethrow(exception);
                    end
                    pause(obj.sleepSeconds(index));
                end
            end
        end

        function value=sleepSeconds(obj,retryIndex)
            value=min(obj.MinDelay*obj.Multiplier^retryIndex,obj.MaxDelay);
        end
    end
end
