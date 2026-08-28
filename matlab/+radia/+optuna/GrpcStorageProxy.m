classdef GrpcStorageProxy < radia.optuna.RDBStorage
    %GRPCSTORAGEPROXY Optuna 4.9 gRPC storage client.

    properties (SetAccess=immutable)
        host (1,1) string
        port (1,1) double
    end

    methods
        function obj=GrpcStorageProxy(options)
            arguments
                options.host (1,1) string = "localhost"
                options.port (1,1) double {mustBeInteger,mustBePositive} = 13000
            end
            module=py.importlib.import_module("optuna.storages");
            proxyClass=py.builtins.getattr(module,"GrpcStorageProxy");
            proxy=proxyClass(pyargs("host",char(options.host), ...
                "port",int64(options.port)));
            obj@radia.optuna.RDBStorage( ...
                "grpc://"+options.host+":"+string(options.port), ...
                InternalPythonStorage=proxy);
            obj.host=options.host;
            obj.port=options.port;
        end

        function close(obj)
            obj.PythonStorage.close();
        end

        function wait_server_ready(obj,timeout)
            if nargin<2 || isempty(timeout)
                obj.PythonStorage.wait_server_ready();
            else
                try
                    obj.PythonStorage.wait_server_ready(double(timeout));
                catch cause
                    if contains(string(cause.message),"GRPC connection timeout")
                        error("radia:optuna:GrpcConnectionTimeout", ...
                            "GRPC connection timeout");
                    end
                    rethrow(cause)
                end
            end
        end
    end
end
