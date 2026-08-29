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
        UseSeparableCMA (1,1) logical = false
        WithMargin (1,1) logical = false
        LrAdapt (1,1) logical = false
        RestartStrategy (1,1) string = ""
        IncPopsize (1,1) double = -1
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
        CandidateRawPoints double = zeros(0,0)
        CandidateRawTrialNumbers double = zeros(0,1)
        SourceTrials cell = cell(0,1)
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
                options.RestartStrategy = []
                options.IncPopsize (1,1) double = -1
                options.UseSeparableCMA (1,1) logical = false
                options.WithMargin (1,1) logical = false
                options.LrAdapt (1,1) logical = false
                options.SourceTrials = []
            end
            if ~isempty(options.RestartStrategy) || options.IncPopsize~=-1
                warning("radia:optuna:FutureWarning", ...
                    "restart_strategy has been deprecated in Optuna 4.4.0, " + ...
                    "falls back to none, and will be removed in 6.0.0.");
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
            if ~isempty(fieldnames(options.X0))
                warning("radia:optuna:FutureWarning", ...
                    "x0 has been deprecated in Optuna 4.9.0 and will be removed in 6.0.0.");
            end
            if isfinite(options.Sigma0)
                warning("radia:optuna:FutureWarning", ...
                    "sigma0 has been deprecated in Optuna 4.9.0 and will be removed in 6.0.0.");
            end
            if ~isempty(options.SourceTrials) && ...
                    (~isempty(fieldnames(options.X0)) || isfinite(sigma))
                error("radia:optuna:CMASourceTrials", ...
                    "SourceTrials cannot be combined with X0, Sigma0, or Sigma.");
            end
            if ~isempty(options.SourceTrials) && options.UseSeparableCMA
                error("radia:optuna:CMASourceTrials", ...
                    "SourceTrials cannot be combined with separable CMA-ES.");
            end
            if options.LrAdapt && ...
                    (options.UseSeparableCMA || options.WithMargin)
                error("radia:optuna:CMALearningRate", ...
                    "LrAdapt cannot be combined with separable CMA-ES or margin.");
            end
            if options.UseSeparableCMA && options.WithMargin
                error("radia:optuna:CMAMargin", ...
                    "Separable CMA-ES and CMA-ES with margin cannot be combined.");
            end
            obj.Seed = radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream = ...
                radia.optuna.internal.NumpyRandomState(obj.Seed);
            if isempty(options.IndependentSampler)
                obj.IndependentSampler = radia.optuna.RandomSampler(obj.Seed);
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
            obj.UseSeparableCMA=options.UseSeparableCMA;
            obj.WithMargin=options.WithMargin;
            obj.LrAdapt=options.LrAdapt;
            if ~isempty(options.RestartStrategy)
                obj.RestartStrategy=string(options.RestartStrategy);
            end
            obj.IncPopsize=options.IncPopsize;
            if isempty(options.SourceTrials)
                obj.SourceTrials=cell(0,1);
            elseif iscell(options.SourceTrials)
                obj.SourceTrials=reshape(options.SourceTrials,[],1);
            else
                obj.SourceTrials=reshape(num2cell(options.SourceTrials),[],1);
            end
            if obj.ConsiderPrunedTrials || obj.UseSeparableCMA || ...
                    obj.WithMargin || obj.LrAdapt || ~isempty(obj.SourceTrials)
                warning("radia:optuna:ExperimentalWarning", ...
                    "The requested advanced CmaEsSampler option is experimental in Optuna 4.9.0.");
            end
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
            if obj.eligibleTrialCount(study) < obj.NStartupTrials
                return
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
            if isempty(searchSpace)
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
            if isa(obj.Engine, ...
                    "radia.optuna.internal.CMAWithMarginEvolutionStrategy")
                [candidate,rawCandidate]=obj.Engine.ask();
                rawRow=size(obj.CandidateRawPoints,1)+1;
                obj.CandidateRawPoints(rawRow,:)=rawCandidate;
                obj.CandidateRawTrialNumbers(rawRow,1)=trial.Number;
            else
                candidate = obj.Engine.ask();
            end
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
            eligible=trial.State=="COMPLETE" || ...
                (obj.ConsiderPrunedTrials && trial.State=="PRUNED" && ...
                ~isempty(trial.IntermediateValues));
            if ~eligible || isempty(obj.Engine) || ...
                    isempty(obj.SearchSpace)
                obj.recordState(study, trial.Number);
                return
            end
            point = zeros(1, numel(obj.SearchSpace));
            rawRow=find(obj.CandidateRawTrialNumbers==trial.Number,1);
            if ~isempty(rawRow)
                if rawRow>size(obj.CandidateRawPoints,1)
                    error("radia:optuna:CMAState", ...
                        "Stored CMA-ES-with-margin candidate state is inconsistent.");
                end
                point=obj.CandidateRawPoints(rawRow,:);
            end
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
                if isempty(rawRow)
                    point(index) = obj.toInternal( ...
                        trial.Params.(key), ...
                        obj.SearchSpace(index).distribution);
                end
            end
            if any(~isfinite(point))
                obj.recordState(study, trial.Number);
                return
            end
            fitness = trial.Value;
            if trial.State=="PRUNED"
                [~,last]=max(trial.IntermediateValues.Step);
                fitness=trial.IntermediateValues.Value(last);
            end
            if study.Directions(1) == "maximize"
                fitness = -fitness;
            end
            obj.PopulationPoints(end+1,:) = point;
            obj.PopulationFitness(end+1,1) = fitness;
            obj.PopulationTrialNumbers(end+1,1) = trial.Number;
            if ~isempty(rawRow)
                obj.CandidateRawPoints(rawRow,:)=[];
                obj.CandidateRawTrialNumbers(rawRow)=[];
            end
            if size(obj.PopulationPoints,1) == obj.Engine.PopulationSize
                obj.Engine.tell(obj.PopulationPoints, obj.PopulationFitness);
                obj.PopulationPoints = zeros(0, numel(obj.SearchSpace));
                obj.PopulationFitness = zeros(0,1);
                obj.PopulationTrialNumbers = zeros(0,1);
                obj.CandidateRawPoints = zeros(0, numel(obj.SearchSpace));
                obj.CandidateRawTrialNumbers = zeros(0,1);
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
                obj.CandidateRawPoints = zeros(0,0);
                obj.CandidateRawTrialNumbers = zeros(0,1);
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
            covariance=zeros(0,0);
            if isempty(obj.SourceTrials)
                if ~isfinite(sigma),sigma=1/6;end
            else
                [mean,sigma,covariance]=obj.warmStart(searchSpace);
            end
            sigma=max(sigma,1e-10);
            engineSeed=randi(obj.Stream,2^31-3);
            bounds=repmat([0,1],dimension,1);
            if obj.UseSeparableCMA && dimension>1
                obj.Engine=radia.optuna.internal. ...
                    SeparableCMAEvolutionStrategy(mean,sigma, ...
                    Bounds=bounds,PopulationSize=obj.PopulationSize, ...
                    Seed=engineSeed,MaxResampling=10*dimension);
            elseif obj.UseSeparableCMA
                warning("radia:optuna:UserWarning", ...
                    "Separable CMA-ES is ignored for a one-dimensional search space.");
                obj.Engine=radia.optuna.internal.CMAEvolutionStrategy( ...
                    mean,sigma,Bounds=bounds,PopulationSize=obj.PopulationSize, ...
                    Seed=engineSeed,Covariance=covariance, ...
                    MaxResampling=10*dimension);
            elseif obj.WithMargin
                steps=zeros(1,dimension);
                for index=1:dimension
                    distribution=searchSpace(index).distribution;
                    if isfinite(distribution.step) && ~distribution.log && ...
                            distribution.low~=distribution.high
                        steps(index)=distribution.step/ ...
                            (distribution.high-distribution.low);
                    elseif isfinite(distribution.step) && ...
                            distribution.low==distribution.high
                        steps(index)=1;
                    end
                end
                obj.Engine=radia.optuna.internal. ...
                    CMAWithMarginEvolutionStrategy(mean,sigma,Bounds=bounds, ...
                    Steps=steps,PopulationSize=obj.PopulationSize, ...
                    Seed=engineSeed,Covariance=covariance, ...
                    MaxResampling=10*dimension);
            else
                obj.Engine=radia.optuna.internal.CMAEvolutionStrategy( ...
                    mean,sigma,Bounds=bounds,PopulationSize=obj.PopulationSize, ...
                    Seed=engineSeed,Covariance=covariance, ...
                    MaxResampling=10*dimension,LrAdapt=obj.LrAdapt);
            end
            obj.SearchSpace = searchSpace;
            obj.SearchSpaceSignature = obj.searchSpaceFingerprint(searchSpace);
            obj.PopulationPoints = zeros(0, dimension);
            obj.PopulationFitness = zeros(0,1);
            obj.PopulationTrialNumbers = zeros(0,1);
            obj.CandidateRawPoints = zeros(0,dimension);
            obj.CandidateRawTrialNumbers = zeros(0,1);
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
                "candidate_raw_points",obj.CandidateRawPoints, ...
                "candidate_raw_trial_numbers",obj.CandidateRawTrialNumbers, ...
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
                schema=string(state.engine.schema);
                switch schema
                    case "radia.optuna.cma-evolution-state.v1"
                        obj.Engine=radia.optuna.internal. ...
                            CMAEvolutionStrategy.fromSnapshot(state.engine);
                    case "radia.optuna.separable-cma-state.v1"
                        obj.Engine=radia.optuna.internal. ...
                            SeparableCMAEvolutionStrategy.fromSnapshot(state.engine);
                    case "radia.optuna.cma-with-margin-state.v1"
                        obj.Engine=radia.optuna.internal. ...
                            CMAWithMarginEvolutionStrategy.fromSnapshot(state.engine);
                    otherwise
                        error("radia:optuna:CMAState", ...
                            "Stored CMA-ES engine type is unsupported.");
                end
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
            if isfield(state,"candidate_raw_points")
                obj.CandidateRawPoints=double(state.candidate_raw_points);
                obj.CandidateRawTrialNumbers=reshape( ...
                    double(state.candidate_raw_trial_numbers),[],1);
            else
                obj.CandidateRawPoints=zeros(0,numel(obj.SearchSpace));
                obj.CandidateRawTrialNumbers=zeros(0,1);
            end
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

        function count=eligibleTrialCount(obj,study)
            trials=study.TrialTable;
            count=sum(trials.State=="COMPLETE");
            if ~obj.ConsiderPrunedTrials,return,end
            pruned=trials.TrialNumber(trials.State=="PRUNED");
            for number=reshape(pruned,1,[])
                if any(study.IntermediateTable.TrialNumber==number)
                    count=count+1;
                end
            end
        end

        function [mean,sigma,covariance]=warmStart(obj,searchSpace)
            points=zeros(0,numel(searchSpace));
            fitness=zeros(0,1);
            for sourceIndex=1:numel(obj.SourceTrials)
                trial=obj.SourceTrials{sourceIndex};
                if ~isa(trial,"radia.optuna.FrozenTrial") || ...
                        ~(trial.State=="COMPLETE" || ...
                        (obj.ConsiderPrunedTrials && trial.State=="PRUNED"))
                    continue
                end
                value=trial.Value;
                if trial.State=="PRUNED"
                    if isempty(trial.IntermediateValues),continue,end
                    [~,last]=max(trial.IntermediateValues.Step);
                    value=trial.IntermediateValues.Value(last);
                end
                if ~isfinite(value),continue,end
                parameterFields=string(fieldnames(trial.Params));
                distributionFields=string(fieldnames(trial.Distributions));
                if numel(parameterFields)~=numel(searchSpace) || ...
                        numel(distributionFields)~=numel(searchSpace)
                    continue
                end
                point=zeros(1,numel(searchSpace));
                compatible=true;
                for index=1:numel(searchSpace)
                    key=matlab.lang.makeValidName(searchSpace(index).name);
                    if ~isfield(trial.Params,key) || ...
                            ~isfield(trial.Distributions,key)
                        compatible=false;
                        break
                    end
                    sourceDistribution=radia.optuna.internal. ...
                        DistributionCodec.normalize(trial.Distributions.(key));
                    if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                            searchSpace(index).distribution,sourceDistribution)
                        compatible=false;
                        break
                    end
                    point(index)=obj.toInternal(trial.Params.(key), ...
                        searchSpace(index).distribution);
                end
                if compatible && all(isfinite(point))
                    points(end+1,:)=point; %#ok<AGROW>
                    fitness(end+1,1)=value; %#ok<AGROW>
                end
            end
            if isempty(points)
                error("radia:optuna:CMASourceTrials", ...
                    "No compatible SourceTrials were supplied.");
            end
            if obj.AttachedStudy.Directions(1)=="maximize"
                fitness=-fitness;
            end
            [~,order]=sort(fitness,"ascend");
            selected=floor(0.1*size(points,1));
            if selected<1
                error("radia:optuna:CMASourceTrials", ...
                    "At least ten compatible SourceTrials are required by warm start.");
            end
            top=points(order(1:selected),:);
            mean=sum(top,1)/selected;
            promising=0.1^2*eye(size(points,2))+ ...
                (top.'*top)/selected-(mean.'*mean);
            determinant=det(promising);
            sigma=determinant^(1/(2*size(points,2)));
            covariance=promising/(determinant^(1/size(points,2)));
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
