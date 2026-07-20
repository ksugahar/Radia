classdef RandomSampler < handle
    %RANDOMSAMPLER Define-by-run random sampler for MATLAB trials.

    properties (SetAccess=private)
        Stream
    end

    methods
        function obj = RandomSampler(seed)
            if nargin < 1 || isempty(seed)
                seed = 0;
            end
            obj.Stream = RandStream("mt19937ar", "Seed", double(seed));
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options)
            %#ok<INUSD>
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
        end

        function value = sampleInteger(obj, study, trial, name, low, high) %#ok<INUSD>
            if ~(isfinite(low) && isfinite(high) && low == floor(low) && ...
                    high == floor(high) && low <= high)
                error("radia:optuna:Bounds", "Integer bounds must be finite integers with low <= high.");
            end
            value = floor(low + rand(obj.Stream, 1, 1) * (high - low + 1));
        end

        function value = sampleCategorical(obj, study, trial, name, choices) %#ok<INUSD>
            if isempty(choices)
                error("radia:optuna:Choices", "Categorical choices must not be empty.");
            end
            if iscell(choices)
                count = numel(choices);
            else
                count = numel(choices);
            end
            index = 1 + floor(rand(obj.Stream, 1, 1) * count);
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
        end

        function beforeTrial(obj, study, trial) %#ok<INUSD>
            % Random sampling has no state update beyond its random stream.
        end

        function afterTrial(obj, study, trial) %#ok<INUSD>
            % Keep the sampler protocol compatible with stateful samplers.
        end
    end
end
