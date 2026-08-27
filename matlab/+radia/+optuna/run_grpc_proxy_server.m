function handle=run_grpc_proxy_server(storage,options)
%RUN_GRPC_PROXY_SERVER Run Optuna's gRPC storage server.
arguments
    storage (1,1) radia.optuna.RDBStorage
    options.host (1,1) string = "localhost"
    options.port (1,1) double {mustBeInteger,mustBePositive} = 13000
    options.thread_pool = []
    options.Background (1,1) logical = false
end
warning("radia:optuna:ExperimentalWarning", ...
    "run_grpc_proxy_server is experimental in Optuna 4.9.0.");
if ~isempty(options.thread_pool)
    error("radia:optuna:GrpcThreadPool", ...
        "thread_pool cannot cross the required MATLAB/Python process boundary.");
end
environment=pyenv;
if environment.Status=="NotLoaded"
    environment=pyenv(ExecutionMode="InProcess");
end
code=[ ...
    "import optuna,sys", ...
    "storage=optuna.storages.RDBStorage(sys.argv[1])", ...
    "optuna.storages.run_grpc_proxy_server(" + ...
    "storage,host=sys.argv[2],port=int(sys.argv[3]))"];
command=java.util.ArrayList();
command.add(char(environment.Executable));
command.add(char("-c"));
command.add(char(strjoin(code,";")));
command.add(char(storage.url));
command.add(char(options.host));
command.add(char(string(options.port)));
builder=java.lang.ProcessBuilder(command);
builder.inheritIO();
process=builder.start();
if options.Background
    handle=radia.optuna.GrpcServerHandle(process);
else
    handle=[];
    exitCode=process.waitFor();
    if exitCode~=0
        error("radia:optuna:GrpcServer", ...
            "Optuna gRPC server exited with status %d.",exitCode);
    end
end
end
