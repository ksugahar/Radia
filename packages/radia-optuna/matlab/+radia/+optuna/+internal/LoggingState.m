classdef LoggingState
    %LOGGINGSTATE Process-local state for the Optuna-compatible logger API.

    methods (Static)
        function value=get()
            value=radia.optuna.internal.LoggingState.access("get",[]);
        end

        function value=set(field,value)
            state=radia.optuna.internal.LoggingState.access("get",[]);
            state.(field)=value;
            value=radia.optuna.internal.LoggingState.access("set",state);
        end
    end

    methods (Static, Access=private)
        function state=access(action,value)
            persistent current
            if isempty(current)
                current=struct("verbosity",20, ...
                    "default_handler",true,"propagation",false);
            end
            if action=="set"
                current=value;
            end
            state=current;
        end
    end
end
