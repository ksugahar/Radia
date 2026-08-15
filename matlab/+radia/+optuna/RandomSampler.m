classdef RandomSampler < handle
    %RANDOMSAMPLER Define-by-run random sampler for MATLAB trials.

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
                seed = 0;
            end
            obj.Seed = double(seed);
            obj.Stream = RandStream("mt19937ar", "Seed", obj.Seed);
        end

        function value = sampleFloat(obj, study, trial, ~, low, high, options)
            if ~isstruct(options) || ~isfield(options, "Log") || ~isfield(options, "Step")
                error("radia:optuna:SamplerOptions", "Sampler options must contain Log and Step fields.");
            end
            if ~(isfinite(low) && isfinite(high) && low < high)
                error("radia:optuna:Bounds", "Float bounds must satisfy low < high.");
            end
            if options.Log && low <= 0
                error("radia:optuna:LogBounds", "Log-uniform bounds must be positive.");
            end
            if isfinite(options.Step) && options.Step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            u = rand(obj.Stream, 1, 1);
            if options.Log
                value = exp(log(low) + u * (log(high) - log(low)));
            else
                value = low + u * (high - low);
            end
            if isfinite(options.Step)
                value = low + round((value - low) / options.Step) * options.Step;
                value = min(max(value, low), high);
            end
            obj.recordState(study, trial.Number);
        end

        function value = sampleInteger(obj, study, trial, ~, low, high)
            if ~(isfinite(low) && isfinite(high) && low == floor(low) && ...
                    high == floor(high) && low <= high)
                error("radia:optuna:Bounds", "Integer bounds must be finite integers with low <= high.");
            end
            value = floor(low + rand(obj.Stream, 1, 1) * (high - low + 1));
            obj.recordState(study, trial.Number);
        end

        function value = sampleCategorical(obj, study, trial, ~, choices)
            if isempty(choices)
                error("radia:optuna:Choices", "Categorical choices must not be empty.");
            end
            count = numel(choices);
            index = 1 + floor(rand(obj.Stream, 1, 1) * count);
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
                obj.Stream = RandStream("mt19937ar", "Seed", obj.Seed);
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
