classdef MOTPESampler < handle
    %MOTPESAMPLER Multi-objective TPE with Optuna-style Parzen mixtures.
    %   Pareto rank and hypervolume contribution select and weight the
    %   below trials. Candidate generation and acquisition evaluation share
    %   the same Parzen implementation as TPESampler.

    properties (SetAccess=private)
        Stream
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.1
        MaxGoodTrials (1,1) double = 25
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
        ConsiderMagicClip (1,1) logical = true
        ConsiderEndpoints (1,1) logical = false
        ConstraintsFcn = []
    end

    methods
        function obj = MOTPESampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.NStartupTrials (1,1) double ...
                    {mustBeInteger, mustBeNonnegative} = 10
                options.Gamma (1,1) double = 0.1
                options.MaxGoodTrials (1,1) double ...
                    {mustBeInteger, mustBePositive} = 25
                options.NumberOfEIChoices (1,1) double ...
                    {mustBeInteger, mustBePositive} = 24
                options.PriorWeight (1,1) double = 1
                options.ConsiderMagicClip (1,1) logical = true
                options.ConsiderEndpoints (1,1) logical = false
                options.ConstraintsFcn = []
            end
            if options.Gamma <= 0 || options.Gamma > 1
                error("radia:optuna:MOTPEGamma", ...
                    "Gamma must be greater than zero and at most one.");
            end
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn, "function_handle")
                error("radia:optuna:ConstraintsFcn", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.Stream = RandStream("mt19937ar", "Seed", options.Seed);
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.MaxGoodTrials = options.MaxGoodTrials;
            obj.NumberOfEIChoices = options.NumberOfEIChoices;
            obj.PriorWeight = options.PriorWeight;
            obj.ConsiderMagicClip = options.ConsiderMagicClip;
            obj.ConsiderEndpoints = options.ConsiderEndpoints;
            obj.ConstraintsFcn = options.ConstraintsFcn;
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options) %#ok<INUSD>
            obj.validate(low, high, options);
            if low == high
                value = low;
                return
            end
            [x, objectives, trialNumbers] = ...
                radia.optuna.internal.ParetoSupport.numericObservations( ...
                study, name);
            valid = isfinite(x) & x >= low & x <= high & ...
                all(isfinite(objectives), 2);
            if options.Log
                valid = valid & x > 0;
            end
            x = x(valid);
            objectives = objectives(valid, :);
            trialNumbers = trialNumbers(valid);
            if numel(x) < obj.NStartupTrials
                value = obj.quantize( ...
                    obj.uniform(low, high, options.Log), ...
                    low, high, options.Step);
                return
            end

            nGood = min(obj.MaxGoodTrials, ceil(obj.Gamma * numel(x)));
            [goodMask, goodWeights] = ...
                radia.optuna.internal.ParetoSupport.splitMOTPE( ...
                study, trialNumbers, objectives, nGood);
            estimatorOptions = { ...
                "Log", options.Log, ...
                "Step", options.Step, ...
                "PriorWeight", obj.PriorWeight, ...
                "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                "ConsiderEndpoints", obj.ConsiderEndpoints};
            below = radia.optuna.internal.ParzenEstimator.numerical( ...
                x(goodMask), low, high, estimatorOptions{:}, ...
                ObservationWeights=goodWeights);
            above = radia.optuna.internal.ParzenEstimator.numerical( ...
                x(~goodMask), low, high, estimatorOptions{:});
            candidates = ...
                radia.optuna.internal.ParzenEstimator.sampleNumerical( ...
                below, obj.Stream, obj.NumberOfEIChoices);
            acquisition = ...
                radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                below, candidates) - ...
                radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                above, candidates);
            [~, best] = max(acquisition);
            value = candidates(best);
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            if low ~= floor(low) || high ~= floor(high) || low > high
                error("radia:optuna:Bounds", ...
                    "Integer bounds must be finite integers with low <= high.");
            end
            if low == high
                value = low;
                return
            end
            value = obj.sampleFloat(study, trial, name, low, high, ...
                struct("Log", false, "Step", 1));
            value = min(max(round(value), low), high);
        end

        function value = sampleCategorical( ...
                obj, study, trial, name, choices) %#ok<INUSD>
            if isempty(choices)
                error("radia:optuna:Choices", ...
                    "Categorical choices must not be empty.");
            end
            [tokens, objectives, trialNumbers] = ...
                radia.optuna.internal.ParetoSupport.categoricalObservations( ...
                study, name);
            choiceTokens = obj.choiceTokens(choices);
            observed = zeros(numel(tokens), 1);
            valid = all(isfinite(objectives), 2);
            for index = 1:numel(tokens)
                match = find(choiceTokens == tokens(index), 1);
                if isempty(match)
                    valid(index) = false;
                else
                    observed(index) = match;
                end
            end
            observed = observed(valid);
            objectives = objectives(valid, :);
            trialNumbers = trialNumbers(valid);
            count = numel(choiceTokens);
            if numel(observed) < obj.NStartupTrials
                index = 1 + floor(rand(obj.Stream, 1, 1) * count);
                value = obj.choiceAt(choices, index);
                return
            end

            nGood = min(obj.MaxGoodTrials, ceil(obj.Gamma * numel(observed)));
            [goodMask, goodWeights] = ...
                radia.optuna.internal.ParetoSupport.splitMOTPE( ...
                study, trialNumbers, objectives, nGood);
            below = radia.optuna.internal.ParzenEstimator.categorical( ...
                observed(goodMask), count, PriorWeight=obj.PriorWeight, ...
                ObservationWeights=goodWeights);
            above = radia.optuna.internal.ParzenEstimator.categorical( ...
                observed(~goodMask), count, PriorWeight=obj.PriorWeight);
            candidates = ...
                radia.optuna.internal.ParzenEstimator.sampleCategorical( ...
                below, obj.Stream, obj.NumberOfEIChoices);
            acquisition = ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                below, candidates) - ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                above, candidates);
            [~, best] = max(acquisition);
            value = obj.choiceAt(choices, candidates(best));
        end

        function beforeTrial(obj, study, trial) %#ok<INUSD>
        end

        function afterTrial(obj, study, trial)
            if trial.State == "COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial, obj.ConstraintsFcn(trial));
            end
        end
    end

    methods (Access=private)
        function validate(~, low, high, options)
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", ...
                    "Float bounds must satisfy low <= high.");
            end
            if options.Log && low <= 0
                error("radia:optuna:LogBounds", ...
                    "Log bounds must be positive.");
            end
            if isfinite(options.Step) && options.Step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            if options.Log && isfinite(options.Step) && low - options.Step / 2 <= 0
                error("radia:optuna:LogBounds", ...
                    "Expanded log-distribution support must be positive.");
            end
        end

        function value = uniform(obj, low, high, logScale)
            u = rand(obj.Stream, 1, 1);
            if logScale
                value = exp(log(low) + u * (log(high) - log(low)));
            else
                value = low + u * (high - low);
            end
        end

        function value = quantize(~, value, low, high, step)
            if isfinite(step)
                value = low + round((value - low) / step) * step;
            end
            value = min(max(value, low), high);
        end

        function tokens = choiceTokens(obj, choices)
            tokens = strings(numel(choices), 1);
            for index = 1:numel(choices)
                tokens(index) = obj.token(obj.choiceAt(choices, index));
            end
        end

        function token = token(~, value)
            if isnumeric(value) && isscalar(value)
                token = "numeric:" + string(value, "%.17g");
            else
                token = string(jsonencode(value));
            end
        end

        function value = choiceAt(~, choices, index)
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
        end
    end
end
