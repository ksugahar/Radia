classdef Environment < handle
    %ENVIRONMENT Small MATLAB-native reset/step environment for Radia solvers.
    %   The callback boundary is deliberately independent of Reinforcement
    %   Learning Toolbox.  Use toFunctionEnv when that toolbox is available.

    properties (SetAccess=private)
        StepFcn
        ResetFcn
        SeedFcn
        MaxSteps (1,1) double = Inf
        Observation = []
        StepCount (1,1) double = 0
        IsDone (1,1) logical = false
        LastInfo struct = struct()
    end

    methods
        function obj = Environment(options)
            arguments
                options.StepFcn (1,1) function_handle
                options.ResetFcn (1,1) function_handle
                options.SeedFcn = []
                options.MaxSteps (1,1) double = Inf
            end
            if ~(isinf(options.MaxSteps) || ...
                    (isfinite(options.MaxSteps) && options.MaxSteps >= 1 && ...
                    options.MaxSteps == floor(options.MaxSteps)))
                error("radia:rl:MaxSteps", ...
                    "MaxSteps must be a positive integer or Inf.");
            end
            if ~isempty(options.SeedFcn) && ~isa(options.SeedFcn, "function_handle")
                error("radia:rl:SeedFcn", "SeedFcn must be a function handle or empty.");
            end
            obj.StepFcn = options.StepFcn;
            obj.ResetFcn = options.ResetFcn;
            obj.SeedFcn = options.SeedFcn;
            obj.MaxSteps = options.MaxSteps;
        end

        function [observation, info] = reset(obj)
            [observation, info] = obj.ResetFcn();
            obj.validateObservation(observation);
            obj.StepCount = 0;
            obj.IsDone = false;
            obj.Observation = observation;
            obj.LastInfo = obj.normalizeInfo(info);
            info = obj.LastInfo;
        end

        function [observation, reward, isDone, info] = step(obj, action)
            if obj.IsDone
                error("radia:rl:EpisodeDone", ...
                    "The episode is complete; call reset before step.");
            end
            [observation, reward, isDone, info] = obj.StepFcn(action);
            obj.validateObservation(observation);
            if ~isnumeric(reward) || ~isscalar(reward) || ~isfinite(reward)
                error("radia:rl:Reward", "StepFcn must return one finite numeric reward.");
            end
            if ~islogical(isDone) || ~isscalar(isDone)
                error("radia:rl:Done", "StepFcn must return one logical isDone value.");
            end
            obj.StepCount = obj.StepCount + 1;
            isDone = logical(isDone || obj.StepCount >= obj.MaxSteps);
            info = obj.normalizeInfo(info);
            if isDone && obj.StepCount >= obj.MaxSteps && ~isfield(info, "termination_reason")
                info.termination_reason = "max_steps";
            end
            obj.Observation = observation;
            obj.IsDone = isDone;
            obj.LastInfo = info;
        end

        function seed(obj, value)
            arguments
                obj
                value (1,1) double {mustBeInteger, mustBeFinite}
            end
            if ~isempty(obj.SeedFcn)
                obj.SeedFcn(value);
            else
                rng(value, "twister");
            end
        end

        function rlEnvironment = toFunctionEnv(obj, observationInfo, actionInfo)
            if exist("rlFunctionEnv", "file") ~= 2 && ...
                    exist("rlFunctionEnv", "builtin") ~= 5
                error("radia:rl:MissingToolbox", ...
                    "Reinforcement Learning Toolbox is required for rlFunctionEnv.");
            end
            rlEnvironment = rlFunctionEnv(observationInfo, actionInfo, ...
                @(action, loggedSignals) obj.stepForRL(action, loggedSignals), ...
                @() obj.resetForRL());
        end
    end

    methods (Access=private)
        function [observation, reward, isDone, loggedSignals] = stepForRL( ...
                obj, action, loggedSignals)
            [observation, reward, isDone, info] = obj.step(action);
            loggedSignals = info;
        end

        function [observation, loggedSignals] = resetForRL(obj)
            [observation, loggedSignals] = obj.reset();
        end

        function validateObservation(~, observation)
            if ~(isnumeric(observation) && ~isempty(observation) && ...
                    all(isfinite(observation), "all"))
                error("radia:rl:Observation", ...
                    "ResetFcn and StepFcn must return finite numeric observations.");
            end
        end

        function info = normalizeInfo(~, info)
            if isempty(info)
                info = struct();
            elseif ~isstruct(info) || ~isscalar(info)
                error("radia:rl:Info", "ResetFcn and StepFcn info must be scalar structs.");
            end
        end
    end
end
