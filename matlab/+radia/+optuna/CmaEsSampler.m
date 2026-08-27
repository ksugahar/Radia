classdef CmaEsSampler < radia.optuna.BaseSampler
    %CMAESSAMPLER Table-persistent full-covariance CMA-ES sampler.
    %   Numeric parameters in the completed-trial intersection search space
    %   are sampled jointly. Categorical and non-intersection parameters use
    %   an independent seeded random proposal.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 1
        PopulationSize (1,1) double = 0
        Sigma0 (1,1) double = NaN
        Sigma (1,1) double = NaN
        X0 struct = struct()
        ConsiderPrunedTrials (1,1) logical = false
        WarnIndependentSampling (1,1) logical = true
    end

    properties (Access=private)
        IndependentSampler
        Engine = []
        SearchSpace struct = ...
            radia.optuna.internal.IntersectionSearchSpace.empty()
        SearchSpaceSignature (1,1) string = ""
        PopulationPoints double = zeros(0,0)
        PopulationFitness double = zeros(0,1)
        PopulationTrialNumbers double = zeros(0,1)
        OptimizerCheckpointed (1,1) logical = false
        AttachedStudy = []
        Restored (1,1) logical = false
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.cma-sampler-state.v1"
        SamplerName = "cmaes"
    end

    methods
        function obj = CmaEsSampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.NStartupTrials (1,1) double = 1
                options.PopulationSize (1,1) double = 0
                options.Sigma0 (1,1) double = NaN
                options.Sigma (1,1) double = NaN
                options.X0 struct = struct()
                options.ConsiderPrunedTrials (1,1) logical = false
                options.IndependentSampler = []
                options.WarnIndependentSampling (1,1) logical = true
            end
            if options.NStartupTrials < 0 || ...
                    options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:CMAStartup", ...
                    "NStartupTrials must be a nonnegative integer.");
            end
            if options.PopulationSize ~= 0 && ...
                    (options.PopulationSize < 2 || ...
                    options.PopulationSize ~= floor(options.PopulationSize))
                error("radia:optuna:CMAPopulation", ...
                    "PopulationSize must be zero or an integer of at least two.");
            end
            if isfinite(options.Sigma0) && isfinite(options.Sigma)
                error("radia:optuna:CMASigma", ...
                    "Specify only Sigma0 or the legacy Sigma option.");
            end
            sigma = options.Sigma0;
            if ~isfinite(sigma)
                sigma = options.Sigma;
            end
            if isfinite(sigma) && sigma <= 0
                error("radia:optuna:CMASigma", ...
                    "Sigma0 must be positive when specified.");
            end
            obj.Seed = radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream = ...
                radia.optuna.internal.NumpyRandomState(obj.Seed);
            if isempty(options.IndependentSampler)
                obj.IndependentSampler = radia.optuna.RandomSampler(options.Seed);
            else
                if ~isobject(options.IndependentSampler)
                    error("radia:optuna:CMAIndependentSampler", ...
                        "IndependentSampler must implement the sampler API.");
                end
                samplerMethods = string(methods(options.IndependentSampler));
                required = ["sampleFloat","sampleInteger", ...
                    "sampleCategorical"];
                if any(~ismember(required,samplerMethods))
                    error("radia:optuna:CMAIndependentSampler", ...
                        "IndependentSampler must implement the sampler API.");
                end
                obj.IndependentSampler = options.IndependentSampler;
            end
            obj.NStartupTrials = options.NStartupTrials;
            obj.PopulationSize = options.PopulationSize;
            obj.Sigma0 = sigma;
            obj.Sigma = sigma;
            obj.X0 = options.X0;
            obj.ConsiderPrunedTrials = options.ConsiderPrunedTrials;
            obj.WarnIndependentSampling = options.WarnIndependentSampling;
        end

        function reseed_rng(obj)
            obj.IndependentSampler.reseed_rng();
        end

        function searchSpace = inferRelativeSearchSpace(obj, study, trial) %#ok<INUSD>
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.calculate( ...
                study, IncludePruned=obj.ConsiderPrunedTrials, ...
                NumericOnly=true);
        end

        function searchSpace = infer_relative_search_space(obj, study, trial)
            if nargin < 3
                trial = [];
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
        end

        function beforeTrial(obj, study, trial)
            if numel(study.Directions) ~= 1
                error("radia:optuna:CMAMultiObjective", ...
                    "CMA-ES supports one objective. Use multivariate TPE for multiple objectives.");
            end
            obj.IndependentSampler.beforeTrial(study,trial);
            obj.attach(study);
            completed = study.TrialTable.State == "COMPLETE";
            if sum(completed) < obj.NStartupTrials
                return
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
            if numel(searchSpace) < 2
                return
            end
            signature = obj.searchSpaceFingerprint(searchSpace);
            initialized=false;
            if isempty(obj.Engine) || signature ~= obj.SearchSpaceSignature
                obj.initializeEngine(searchSpace);
                initialized=true;
            elseif ~obj.OptimizerCheckpointed
                % Until the first complete population is told, upstream
                % cannot restore a serialized cmaes optimizer and creates a
                % fresh one on every trial.  Its constructor seed is later
                % overwritten, but the sampler-level randint is observable.
                randi(obj.Stream,2^31-3);
            end
            % Optuna re-seeds the cmaes optimizer for every trial from a
            % separate sampler-level RandomState, then adds trial.number.
            candidateSeed=randi(obj.Stream,2^16-1)+trial.Number;
            obj.Engine.reseed(candidateSeed);
            candidate = obj.Engine.ask();
            if obj.Engine.Generation>0
                obj.OptimizerCheckpointed=true;
            elseif initialized
                obj.OptimizerCheckpointed=false;
            end
            values = cell(1, numel(searchSpace));
            for index = 1:numel(searchSpace)
                values{index} = obj.fromInternal( ...
                    candidate(index), searchSpace(index).distribution);
            end
            trial.setRelativeParameters(searchSpace, values, "cmaes");
            obj.recordState(study, trial.Number);
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options)
            obj.warnIndependent(study, trial, name);
            value=obj.IndependentSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            obj.warnIndependent(study, trial, name);
            value=obj.IndependentSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value = sampleCategorical(obj, study, trial, name, choices)
            obj.warnIndependent(study, trial, name);
            value=obj.IndependentSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function afterTrial(obj, study, trial)
            obj.IndependentSampler.afterTrial(study,trial);
            obj.attach(study);
            if trial.State ~= "COMPLETE" || isempty(obj.Engine) || ...
                    isempty(obj.SearchSpace)
                obj.recordState(study, trial.Number);
                return
            end
            point = zeros(1, numel(obj.SearchSpace));
            for index = 1:numel(obj.SearchSpace)
                key = matlab.lang.makeValidName(obj.SearchSpace(index).name);
                if ~isfield(trial.Params, key)
                    obj.recordState(study, trial.Number);
                    return
                end
                parameterRow = study.ParamTable.TrialNumber == trial.Number & ...
                    study.ParamTable.Name == obj.SearchSpace(index).name;
                if sum(parameterRow) ~= 1
                    obj.recordState(study, trial.Number);
                    return
                end
                stored = radia.optuna.internal.DistributionCodec.decode( ...
                    study.ParamTable.Kind(parameterRow), ...
                    study.ParamTable.Distribution(parameterRow));
                if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                        obj.SearchSpace(index).distribution, stored)
                    obj.recordState(study, trial.Number);
                    return
                end
                point(index) = obj.toInternal( ...
                    trial.Params.(key), ...
                    obj.SearchSpace(index).distribution);
            end
            if any(~isfinite(point))
                obj.recordState(study, trial.Number);
                return
            end
            fitness = trial.Value;
            if study.Directions(1) == "maximize"
                fitness = -fitness;
            end
            obj.PopulationPoints(end+1,:) = point;
            obj.PopulationFitness(end+1,1) = fitness;
            obj.PopulationTrialNumbers(end+1,1) = trial.Number;
            if size(obj.PopulationPoints,1) == obj.Engine.PopulationSize
                obj.Engine.tell(obj.PopulationPoints, obj.PopulationFitness);
                obj.PopulationPoints = zeros(0, numel(obj.SearchSpace));
                obj.PopulationFitness = zeros(0,1);
                obj.PopulationTrialNumbers = zeros(0,1);
            end
            obj.Sigma = obj.Engine.Sigma;
            obj.recordState(study, trial.Number);
        end
    end

    methods (Access=private)
        function warnIndependent(obj, study, trial, name)
            if ~obj.WarnIndependentSampling
                return
            end
            count = sum(study.TrialTable.State == "COMPLETE");
            if obj.ConsiderPrunedTrials
                prunedNumbers = study.TrialTable.TrialNumber( ...
                    study.TrialTable.State == "PRUNED");
                for number = reshape(prunedNumbers,1,[])
                    if any(study.IntermediateTable.TrialNumber == number)
                        count = count + 1;
                    end
                end
            end
            if count < obj.NStartupTrials
                return
            end
            warning("radia:optuna:CMAIndependentSampling", ...
                "Parameter '%s' in trial %d is sampled independently by " + ...
                "%s because dynamic spaces and categorical distributions " + ...
                "are outside CmaEsSampler's relative search space.", ...
                string(name), trial.Number, class(obj.IndependentSampler));
        end

        function attach(obj, study)
            changed = isempty(obj.AttachedStudy) || ...
                ~isequal(obj.AttachedStudy, study);
            if changed
                obj.AttachedStudy = study;
                obj.Stream = ...
                    radia.optuna.internal.NumpyRandomState(obj.Seed);
                obj.Engine = [];
                obj.SearchSpace = ...
                    radia.optuna.internal.IntersectionSearchSpace.empty();
                obj.SearchSpaceSignature = "";
                obj.PopulationPoints = zeros(0,0);
                obj.PopulationFitness = zeros(0,1);
                obj.PopulationTrialNumbers = zeros(0,1);
                obj.OptimizerCheckpointed = false;
                obj.Restored = false;
            end
            if obj.Restored
                return
            end
            state = study.samplerState(obj.SamplerName, obj.StateSchema);
            if ~isempty(state)
                obj.restoreState(state);
            end
            obj.Restored = true;
        end

        function initializeEngine(obj, searchSpace)
            dimension = numel(searchSpace);
            mean = 0.5 * ones(1, dimension);
            fields = fieldnames(obj.X0);
            if ~isempty(fields)
                for index = 1:dimension
                    key = matlab.lang.makeValidName(searchSpace(index).name);
                    if isfield(obj.X0, key)
                        mean(index) = obj.toInternal( ...
                            obj.X0.(key), searchSpace(index).distribution);
                    end
                end
            end
            if any(~isfinite(mean) | mean < 0 | mean > 1)
                error("radia:optuna:CMAX0", ...
                    "X0 values must lie inside their parameter distributions.");
            end
            sigma = obj.Sigma0;
            if ~isfinite(sigma)
                sigma = 1 / 6;
            end
            engineSeed=randi(obj.Stream,2^31-3);
            obj.Engine = ...
                radia.optuna.internal.CMAEvolutionStrategy(mean, sigma, ...
                Bounds=repmat([0,1], dimension, 1), ...
                PopulationSize=obj.PopulationSize, Seed=engineSeed, ...
                MaxResampling=10*dimension);
            obj.SearchSpace = searchSpace;
            obj.SearchSpaceSignature = obj.searchSpaceFingerprint(searchSpace);
            obj.PopulationPoints = zeros(0, dimension);
            obj.PopulationFitness = zeros(0,1);
            obj.PopulationTrialNumbers = zeros(0,1);
            obj.OptimizerCheckpointed = false;
            obj.Sigma = obj.Engine.Sigma;
        end

        function state = snapshot(obj)
            engineState = [];
            generation = 0;
            if ~isempty(obj.Engine)
                engineState = obj.Engine.snapshot();
                generation = obj.Engine.Generation;
            end
            state = struct( ...
                "schema", obj.StateSchema, ...
                "seed", obj.Seed, ...
                "search_space", obj.SearchSpace, ...
                "search_space_signature", obj.SearchSpaceSignature, ...
                "engine", engineState, ...
                "population_points", obj.PopulationPoints, ...
                "population_fitness", obj.PopulationFitness, ...
                "population_trial_numbers", obj.PopulationTrialNumbers, ...
                "optimizer_checkpointed",obj.OptimizerCheckpointed, ...
                "independent_random_state", obj.Stream.State, ...
                "generation", generation);
        end

        function restoreState(obj, state)
            required = ["schema","search_space","search_space_signature", ...
                "engine","population_points","population_fitness", ...
                "population_trial_numbers","independent_random_state"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state, required)) || ...
                    string(state.schema) ~= obj.StateSchema
                error("radia:optuna:CMAState", ...
                    "Stored CMA-ES sampler state is invalid or unsupported.");
            end
            obj.SearchSpace = state.search_space;
            obj.SearchSpaceSignature = string(state.search_space_signature);
            if isempty(state.engine)
                obj.Engine = [];
            else
                obj.Engine = ...
                    radia.optuna.internal.CMAEvolutionStrategy.fromSnapshot( ...
                    state.engine);
                if obj.PopulationSize ~= 0 && ...
                        obj.Engine.PopulationSize ~= obj.PopulationSize
                    error("radia:optuna:CMAState", ...
                        "Stored CMA-ES population size does not match the sampler.");
                end
                obj.Sigma = obj.Engine.Sigma;
            end
            obj.PopulationPoints = double(state.population_points);
            obj.PopulationFitness = reshape( ...
                double(state.population_fitness), [], 1);
            obj.PopulationTrialNumbers = reshape( ...
                double(state.population_trial_numbers), [], 1);
            if isfield(state,"optimizer_checkpointed")
                obj.OptimizerCheckpointed=logical(state.optimizer_checkpointed);
            else
                obj.OptimizerCheckpointed=~isempty(obj.Engine) && ...
                    obj.Engine.Generation>0 && ...
                    ~isempty(obj.PopulationTrialNumbers);
            end
            obj.Stream.State = state.independent_random_state;
        end

        function recordState(obj, study, trialNumber)
            generation = 0;
            if ~isempty(obj.Engine)
                generation = obj.Engine.Generation;
            end
            study.recordSamplerState(obj.SamplerName, obj.StateSchema, ...
                trialNumber, generation, obj.snapshot());
        end

        function signature = searchSpaceFingerprint(~, searchSpace)
            parts = strings(1, numel(searchSpace));
            for index = 1:numel(searchSpace)
                parts(index) = searchSpace(index).name + "=" + ...
                    radia.optuna.internal.DistributionCodec.encode( ...
                    searchSpace(index).distribution);
            end
            signature = strjoin(parts, "|");
        end

        function value = toInternal(~, value, distribution)
            value = double(value);
            low = distribution.low;
            high = distribution.high;
            if distribution.log
                if value <= 0
                    value = NaN;
                    return
                end
                value = log(value);
                low = log(low);
                high = log(high);
            end
            value = (value - low) / (high - low);
            if value < -1e-12 || value > 1 + 1e-12
                value = NaN;
            else
                value = min(max(value, 0), 1);
            end
        end

        function value = fromInternal(obj, value, distribution)
            value = min(max(double(value), 0), 1);
            low = distribution.low;
            high = distribution.high;
            if distribution.log
                value = exp(log(low) + value * (log(high) - log(low)));
            else
                value = low + value * (high - low);
            end
            value = obj.quantize( ...
                value, distribution.low, distribution.high, ...
                distribution.step);
            if distribution.kind == "integer"
                value = round(value);
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

        function validateBounds(~, low, high, logScale, step)
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", ...
                    "Float bounds must satisfy low <= high.");
            end
            if logScale && low <= 0
                error("radia:optuna:LogBounds", ...
                    "Log-uniform bounds must be positive.");
            end
            if isfinite(step) && step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
        end
    end
end
