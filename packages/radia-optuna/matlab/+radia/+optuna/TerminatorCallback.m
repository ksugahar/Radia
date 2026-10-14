classdef TerminatorCallback
    %TERMINATORCALLBACK Study.optimize callback wrapping Terminator.

    properties (SetAccess=private)
        Terminator
    end

    methods
        function obj=TerminatorCallback(terminator)
            arguments
                terminator = radia.optuna.Terminator()
            end
            obj.Terminator=terminator;
        end

        function invoke(obj,study,~)
            if obj.Terminator.shouldTerminate(study), study.stop(); end
        end

        function callback=callback(obj)
            callback=@obj.invoke;
        end
    end
end
