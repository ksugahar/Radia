classdef RandomSampler < radia.optuna.BaseSampler
    %RANDOMSAMPLER Optuna 4.9 RandomSampler for MATLAB trials.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.random-sampler-state.v1"
        SamplerName = "random"
    end

    methods
        function obj = RandomSampler(seed)
            if nargin < 1 || isempty(seed)
                seed = double.empty(1,0);
            end
            obj.Seed = radia.optuna.internal.resolveSeed(seed);
            obj.Stream = radia.optuna.internal.NumpyRandomState(obj.Seed);
        end

        function value = sampleFloat(obj, study, trial, ~, low, high, options)
            if ~isstruct(options) || ~isfield(options, "Log") || ~isfield(options, "Step")
                error("radia:optuna:SamplerOptions", "Sampler options must contain Log and Step fields.");
            end
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", "Float bounds must satisfy low <= high.");
            end
            if options.Log && low <= 0
                error("radia:optuna:LogBounds", "Log-uniform bounds must be positive.");
            end
            if isfinite(options.Step) && options.Step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            if low == high
                value = low;
                return
            end
            u = rand(obj.Stream, 1, 1);
            if options.Log && isfinite(options.Step)
                % This path represents IntDistribution(log=true), which is
                % routed through sampleFloat by Trial.suggest_int.
                transformedLow = log(low - 0.5 * options.Step);
                transformedHigh = log(high + 0.5 * options.Step);
                proposal = exp(transformedLow + u * ...
                    (transformedHigh - transformedLow));
                value = round(proposal);
            elseif options.Log
                value = exp(log(low) + u * (log(high) - log(low)));
            elseif isfinite(options.Step)
                transformedLow = low - 0.5 * options.Step;
                transformedHigh = high + 0.5 * options.Step;
                proposal = transformedLow + u * ...
                    (transformedHigh - transformedLow);
                value = low + round((proposal - low) / options.Step) * ...
                    options.Step;
            else
                value = low + u * (high - low);
            end
            value = min(max(value, low), high);
            obj.recordState(study, trial.Number);
        end

        function value = sampleInteger(obj, study, trial, ~, low, high)
            if ~(isfinite(low) && isfinite(high) && low == floor(low) && ...
                    high == floor(high) && low <= high)
                error("radia:optuna:Bounds", "Integer bounds must be finite integers with low <= high.");
            end
            if low == high
                value = low;
                return
            end
            proposal = (low - 0.5) + rand(obj.Stream, 1, 1) * ...
                (high - low + 1);
            value = min(max(round(proposal), low), high);
            obj.recordState(study, trial.Number);
        end

        function value = sampleCategorical(obj, study, trial, ~, choices)
            if isempty(choices)
                error("radia:optuna:Choices", "Categorical choices must not be empty.");
            end
            count = numel(choices);
            % _SearchSpaceTransform uses one independently sampled one-hot
            % score per choice and selects the first maximum.
            [~, index] = max(rand(obj.Stream, count, 1));
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
            obj.recordState(study, trial.Number);
        end

        function beforeTrial(obj, study, ~)
            obj.attach(study);
        end

        function afterTrial(obj, study, trial)
            obj.recordState(study, trial.Number);
        end
    end

    methods (Access=private)
        function attach(obj, study)
            changed = isempty(obj.AttachedStudy) || ...
                ~isequal(obj.AttachedStudy, study);
            if changed
                obj.AttachedStudy = study;
                obj.Stream = ...
                    radia.optuna.internal.NumpyRandomState(obj.Seed);
                obj.Restored = false;
            end
            if obj.Restored
                return
            end
            state = study.samplerState(obj.SamplerName, obj.StateSchema);
            if ~isempty(state)
                if ~isstruct(state) || ~isscalar(state) || ...
                        ~isfield(state, "schema") || ...
                        ~isfield(state, "seed") || ...
                        ~isfield(state, "random_state") || ...
                        string(state.schema) ~= obj.StateSchema || ...
                        double(state.seed) ~= obj.Seed
                    error("radia:optuna:RandomState", ...
                        "Stored random-sampler state is invalid or incompatible.");
                end
                obj.Stream.State = state.random_state;
            end
            obj.Restored = true;
        end

        function recordState(obj, study, trialNumber)
            if strlength(study.StoragePath)==0
                return
            end
            obj.attach(study);
            state = struct( ...
                "schema", obj.StateSchema, ...
                "seed", obj.Seed, ...
                "random_state", obj.Stream.State);
            study.recordSamplerState(obj.SamplerName, obj.StateSchema, ...
                trialNumber, 0, state);
        end
    end
end
