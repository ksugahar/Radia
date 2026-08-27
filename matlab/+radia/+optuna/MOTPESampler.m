classdef MOTPESampler < radia.optuna.BaseSampler
    %MOTPESAMPLER Multi-objective TPE with Optuna-style Parzen mixtures.
    %   Pareto rank and hypervolume contribution select and weight the
    %   below trials. Candidate generation and acquisition evaluation share
    %   the same Parzen implementation as TPESampler.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.1
        MaxGoodTrials (1,1) double = 25
        GammaFcn = []
        WeightsFcn = []
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
        ConsiderMagicClip (1,1) logical = true
        ConsiderEndpoints (1,1) logical = false
        ConstraintsFcn = []
        CategoricalDistanceFcn = []
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
        IndependentSampler
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.motpe-sampler-state.v1"
        SamplerName = "motpe"
    end

    methods
        function obj = MOTPESampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.NStartupTrials (1,1) double ...
                    {mustBeInteger, mustBeNonnegative} = 10
                options.Gamma (1,1) double = 0.1
                options.MaxGoodTrials (1,1) double ...
                    {mustBeInteger, mustBePositive} = 25
                options.GammaFcn = []
                options.WeightsFcn = []
                options.NumberOfEIChoices (1,1) double ...
                    {mustBeInteger, mustBePositive} = 24
                options.PriorWeight (1,1) double = 1
                options.ConsiderMagicClip (1,1) logical = true
                options.ConsiderEndpoints (1,1) logical = false
                options.ConstraintsFcn = []
                options.CategoricalDistanceFcn = []
            end
            if options.Gamma <= 0 || options.Gamma > 1
                error("radia:optuna:MOTPEGamma", ...
                    "Gamma must be greater than zero and at most one.");
            end
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            if ~isempty(options.GammaFcn) && ...
                    ~isa(options.GammaFcn, "function_handle")
                error("radia:optuna:TPEGamma", ...
                    "GammaFcn must be a function handle.");
            end
            if ~isempty(options.WeightsFcn) && ...
                    ~isa(options.WeightsFcn, "function_handle")
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn must be a function handle.");
            end
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn, "function_handle")
                error("radia:optuna:ConstraintsFcn", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.Seed = radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream = ...
                radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.IndependentSampler = radia.optuna.RandomSampler(options.Seed);
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.MaxGoodTrials = options.MaxGoodTrials;
            obj.GammaFcn = options.GammaFcn;
            obj.WeightsFcn = options.WeightsFcn;
            obj.NumberOfEIChoices = options.NumberOfEIChoices;
            obj.PriorWeight = options.PriorWeight;
            obj.ConsiderMagicClip = options.ConsiderMagicClip;
            obj.ConsiderEndpoints = options.ConsiderEndpoints;
            obj.ConstraintsFcn = options.ConstraintsFcn;
            radia.optuna.internal.CategoricalDistance.validate( ...
                options.CategoricalDistanceFcn);
            obj.CategoricalDistanceFcn=options.CategoricalDistanceFcn;
        end

        function reseed_rng(obj)
            freshSeed=radia.optuna.internal.resolveSeed([]);
            obj.Stream=radia.optuna.internal.NumpyRandomState(freshSeed);
            obj.IndependentSampler.reseed_rng();
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options)
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
                value = obj.IndependentSampler.sampleFloat( ...
                    study, trial, name, low, high, options);
                return
            end

            nGood = obj.goodTrialCount(numel(x));
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
                x(~goodMask), low, high, estimatorOptions{:}, ...
                ObservationWeights=obj.observationWeights(sum(~goodMask)));
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
            obj.recordState(study, trial.Number);
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
                obj, study, trial, name, choices)
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
                value = obj.IndependentSampler.sampleCategorical( ...
                    study, trial, name, choices);
                return
            end

            nGood = obj.goodTrialCount(numel(observed));
            [goodMask, goodWeights] = ...
                radia.optuna.internal.ParetoSupport.splitMOTPE( ...
                study, trialNumbers, objectives, nGood);
            distanceFcn=radia.optuna.internal.CategoricalDistance. ...
                get(obj.CategoricalDistanceFcn,name);
            below = radia.optuna.internal.ParzenEstimator.categorical( ...
                observed(goodMask), count, PriorWeight=obj.PriorWeight, ...
                ObservationWeights=goodWeights,DistanceFcn=distanceFcn, ...
                Choices=choices);
            above = radia.optuna.internal.ParzenEstimator.categorical( ...
                observed(~goodMask), count, PriorWeight=obj.PriorWeight, ...
                ObservationWeights=obj.observationWeights(sum(~goodMask)), ...
                DistanceFcn=distanceFcn,Choices=choices);
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
            obj.recordState(study, trial.Number);
        end

        function beforeTrial(obj, study, trial)
            obj.IndependentSampler.beforeTrial(study,trial);
            obj.attach(study);
        end

        function afterTrial(obj, study, trial)
            if trial.State == "COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial, obj.ConstraintsFcn(trial));
            end
            obj.recordState(study, trial.Number);
        end
    end

    methods (Access=private)
        function count = goodTrialCount(obj, finishedCount)
            if isempty(obj.GammaFcn)
                count = min(obj.MaxGoodTrials, ...
                    ceil(obj.Gamma * finishedCount));
                count = max(1, min(finishedCount, count));
                return
            end
            count = obj.GammaFcn(finishedCount);
            if ~(isnumeric(count) && isreal(count) && isscalar(count) && ...
                    isfinite(count) && count == floor(count) && ...
                    count >= 0 && count <= finishedCount)
                error("radia:optuna:TPEGamma", ...
                    "GammaFcn(%d) must return an integer from 0 through %d.", ...
                    finishedCount, finishedCount);
            end
            count = double(count);
        end

        function weights = observationWeights(obj, count)
            if isempty(obj.WeightsFcn)
                weights = zeros(0,1);
                return
            end
            weights = obj.WeightsFcn(count);
            if ~(isnumeric(weights) && isreal(weights))
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) must return a real numeric vector.", count);
            end
            weights = reshape(double(weights),[],1);
            if numel(weights) ~= count
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) returned %d weights; expected %d.", ...
                    count, numel(weights), count);
            end
            if any(~isfinite(weights)) || any(weights < 0) || ...
                    (count > 0 && sum(weights) <= 0)
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) must return finite nonnegative weights " + ...
                    "with positive total mass.", count);
            end
        end

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
                    error("radia:optuna:MOTPEState", ...
                        "Stored MOTPE sampler state is invalid or incompatible.");
                end
                obj.Stream.State = state.random_state;
            end
            obj.Restored = true;
        end

        function recordState(obj, study, trialNumber)
            obj.attach(study);
            if strlength(study.StoragePath)==0
                return
            end
            state = struct( ...
                "schema", obj.StateSchema, ...
                "seed", obj.Seed, ...
                "random_state", obj.Stream.State);
            study.recordSamplerState(obj.SamplerName, obj.StateSchema, ...
                trialNumber, 0, state);
        end

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
