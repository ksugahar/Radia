classdef JournalRedisBackend < radia.optuna.BaseJournalBackend
    %JOURNALREDISBACKEND Redis-backed Optuna 4.9 journal log.

    properties (SetAccess=immutable)
        url (1,1) string
        use_cluster (1,1) logical
        prefix (1,1) string
    end

    properties (Access=private)
        Client
        PythonBackend
    end

    methods
        function obj=JournalRedisBackend(url,options)
            arguments
                url (1,1) string
                options.use_cluster (1,1) logical = false
                options.prefix (1,1) string = ""
                options.Client = []
            end
            obj.url=url;
            obj.use_cluster=options.use_cluster;
            obj.prefix=options.prefix;
            obj.Client=options.Client;
            if isempty(obj.Client)
                try
                    module=py.importlib.import_module( ...
                        "optuna.storages.journal");
                    backendClass=py.builtins.getattr( ...
                        module,"JournalRedisBackend");
                    obj.PythonBackend=backendClass(char(url),pyargs( ...
                        "use_cluster",options.use_cluster, ...
                        "prefix",char(options.prefix)));
                catch cause
                    error("radia:optuna:JournalRedisPython", ...
                        "JournalRedisBackend requires the Optuna redis extra: %s", ...
                        cause.message);
                end
            end
        end

        function append_logs(obj,logs)
            entries=obj.logCells(logs);
            if ~isempty(obj.PythonBackend)
                jsonModule=py.importlib.import_module("json");
                loads=py.builtins.getattr(jsonModule,"loads");
                pythonLogs=cellfun(@(entry)loads(char(jsonencode(entry))), ...
                    entries,"UniformOutput",false);
                obj.PythonBackend.append_logs(py.list(pythonLogs));
                return
            end
            numberKey=obj.prefix+":log_number";
            obj.invoke("setnx",numberKey,-1);
            for index=1:numel(entries)
                number=double(obj.invoke("incr",numberKey,1));
                obj.invoke("set",obj.logKey(number),jsonencode(entries{index}));
            end
        end

        function logs=read_logs(obj,log_number_from)
            arguments
                obj
                log_number_from (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            if ~isempty(obj.PythonBackend)
                source=cell(py.list( ...
                    obj.PythonBackend.read_logs(int64(log_number_from))));
                logs=cell(numel(source),1);
                jsonModule=py.importlib.import_module("json");
                dumps=py.builtins.getattr(jsonModule,"dumps");
                for index=1:numel(source)
                    logs{index}=jsondecode(char(string(dumps(source{index}))));
                end
                return
            end
            maximum=obj.invoke("get",obj.prefix+":log_number");
            if isempty(maximum)
                logs=cell(0,1);
                return
            end
            maximum=str2double(string(maximum));
            logs=cell(max(0,maximum-log_number_from+1),1);
            destination=0;
            for number=log_number_from:maximum
                encoded=obj.invoke("get",obj.logKey(number));
                if isempty(encoded)
                    error("radia:optuna:JournalRedisMissing", ...
                        "Redis journal log %d is missing.",number);
                end
                destination=destination+1;
                logs{destination}=jsondecode(char(string(encoded)));
            end
        end

        function save_snapshot(obj,snapshot)
            bytes=reshape(uint8(snapshot),1,[]);
            if ~isempty(obj.PythonBackend)
                obj.PythonBackend.save_snapshot( ...
                    py.bytes(py.list(num2cell(bytes))));
            else
                obj.invoke("set",obj.prefix+":snapshot",bytes);
            end
        end

        function snapshot=load_snapshot(obj)
            if ~isempty(obj.PythonBackend)
                source=obj.PythonBackend.load_snapshot();
                if isa(source,"py.NoneType")
                    snapshot=[];
                else
                    snapshot=uint8(source);
                end
            else
                snapshot=obj.invoke("get",obj.prefix+":snapshot");
                if ~isempty(snapshot), snapshot=reshape(uint8(snapshot),1,[]); end
            end
        end

        function backend=pythonBackend(obj)
            if isempty(obj.PythonBackend)
                error("radia:optuna:JournalRedisPython", ...
                    "A MATLAB-injected Redis client cannot be passed to " + ...
                    "Optuna JournalStorage replay.");
            end
            backend=obj.PythonBackend;
        end
    end

    methods (Access=private)
        function value=invoke(obj,name,varargin)
            if ~isstruct(obj.Client) || ~isfield(obj.Client,name) || ...
                    ~isa(obj.Client.(name),"function_handle")
                error("radia:optuna:JournalRedisClient", ...
                    "Redis client does not implement '%s'.",name);
            end
            value=obj.Client.(name)(varargin{:});
        end

        function key=logKey(obj,number)
            key=obj.prefix+":log:"+string(number);
        end

        function entries=logCells(~,logs)
            if iscell(logs)
                entries=reshape(logs,[],1);
            elseif isstruct(logs)
                entries=num2cell(reshape(logs,[],1));
            else
                error("radia:optuna:JournalLogs", ...
                    "logs must be a struct array or cell array of structs.");
            end
        end
    end
end
