classdef GrpcServerHandle < handle
    %GRPCSERVERHANDLE Lifecycle owner for a background Optuna gRPC server.

    properties (Access=private)
        Process
        Stopped (1,1) logical = false
    end

    methods (Hidden=true)
        function obj=GrpcServerHandle(process)
            obj.Process=process;
        end
    end

    methods
        function stop(obj,grace) %#ok<INUSD>
            if obj.Stopped, return, end
            if obj.Process.isAlive()
                obj.Process.destroy();
                pause(0.1);
            end
            if obj.Process.isAlive()
                obj.Process.destroyForcibly();
            end
            obj.Stopped=true;
        end

        function delete(obj)
            if ~obj.Stopped
                try
                    obj.stop(0);
                catch
                end
            end
        end
    end
end
