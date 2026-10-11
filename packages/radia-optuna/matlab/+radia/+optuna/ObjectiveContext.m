classdef ObjectiveContext < handle
    %OBJECTIVECONTEXT Trial services exposed to a model-shaped objective.
    %
    % A callback used by OptimizationSession may accept (values, context).
    % The values struct remains independent of Optuna. Context supplies the
    % intermediate-value and pruning operations that need the live Trial.

    properties (SetAccess=private)
        Trial = []
    end

    properties
        UserData = []
    end

    methods
        function obj = ObjectiveContext(trial)
            arguments
                trial (1,1) radia.optuna.Trial
            end
            obj.Trial = trial;
        end

        function report(obj, value, step)
            obj.Trial.report(value, step);
        end

        function decision = shouldPrune(obj)
            decision = obj.Trial.should_prune();
        end

        function decision = should_prune(obj)
            decision = obj.shouldPrune();
        end

        function prune(obj)
            error("radia:optuna:TrialPruned", ...
                "Trial %d was pruned by the objective context.", ...
                obj.Trial.Number);
        end

        function setUserAttr(obj, name, value)
            obj.Trial.set_user_attr(name, value);
        end

        function set_user_attr(obj, name, value)
            obj.setUserAttr(name, value);
        end
    end
end
