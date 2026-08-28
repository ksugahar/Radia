classdef NSGAIISampler < radia.optuna.BaseGASampler
    %NSGAIISAMPLER Optuna-compatible generational constrained NSGA-II.
    %   The default is UniformCrossover. All Optuna 4.9 built-in numerical
    %   crossovers are available under radia.optuna.nsgaii. Categorical
    %   parameters always use uniform crossover, and mutated/dynamic
    %   parameters use the seeded independent random fallback.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
        PopulationSize (1,1) double = 50
        MutationProbability (1,1) double = NaN
        CrossoverProbability (1,1) double = 0.9
        SwappingProbability (1,1) double = 0.5
        Crossover = []
        % Retained for source compatibility; random fallback mutation does
        % not use the former Gaussian MutationScale.
        MutationScale (1,1) double = 0.1
        ConstraintsFcn = []
        ElitePopulationSelectionStrategy = []
        ChildGenerationStrategy = []
        AfterTrialStrategy = []
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
        IndependentSampler
        GenerationTrialNumbers double = zeros(0,1)
        GenerationAssignments double = zeros(0,1)
        ParentCaches struct = struct( ...
            "generation",{},"trial_numbers",{})
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.nsgaii-sampler-state.v3"
        LegacyStateSchemaV2 = "radia.optuna.nsgaii-sampler-state.v2"
        LegacyStateSchemaV1 = "radia.optuna.nsgaii-sampler-state.v1"
        SamplerName = "nsgaii"
    end

    methods
        function obj = NSGAIISampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.PopulationSize (1,1) double ...
                    {mustBeInteger,mustBePositive} = 50
                options.MutationProbability (1,1) double = NaN
                options.Crossover = []
                options.CrossoverProbability (1,1) double = 0.9
                options.SwappingProbability (1,1) double = 0.5
                options.MutationScale (1,1) double {mustBePositive} = 0.1
                options.ConstraintsFcn = []
                options.ElitePopulationSelectionStrategy = []
                options.ChildGenerationStrategy = []
                options.AfterTrialStrategy = []
            end
            obj@radia.optuna.BaseGASampler(options.PopulationSize);
            mutationProbability = options.MutationProbability;
            if ~(isnan(mutationProbability) || ...
                    (isfinite(mutationProbability) && ...
                    mutationProbability >= 0 && mutationProbability <= 1)) || ...
                    ~isfinite(options.CrossoverProbability) || ...
                    options.CrossoverProbability < 0 || ...
                    options.CrossoverProbability > 1 || ...
                    ~isfinite(options.SwappingProbability) || ...
                    options.SwappingProbability < 0 || ...
                    options.SwappingProbability > 1
                error("radia:optuna:NSGAIIProbability", ...
                    "Mutation probability must be NaN or in [0,1], and " + ...
                    "crossover/swapping probabilities must be in [0,1].");
            end
            functionOptions = {options.ConstraintsFcn, ...
                options.ElitePopulationSelectionStrategy, ...
                options.ChildGenerationStrategy,options.AfterTrialStrategy};
            functionNames = ["ConstraintsFcn", ...
                "ElitePopulationSelectionStrategy", ...
                "ChildGenerationStrategy","AfterTrialStrategy"];
            for index = 1:numel(functionOptions)
                if ~isempty(functionOptions{index}) && ...
                        ~isa(functionOptions{index},"function_handle")
                    error("radia:optuna:NSGAIIStrategy", ...
                        "%s must be a function handle.",functionNames(index));
                end
            end
            crossover = options.Crossover;
            if isempty(crossover)
                crossover = radia.optuna.nsgaii.UniformCrossover( ...
                    SwappingProbability=options.SwappingProbability);
            end
            if ~isa(crossover,"radia.optuna.nsgaii.BaseCrossover")
                error("radia:optuna:NSGAIICrossover", ...
                    "Crossover must derive from radia.optuna.nsgaii.BaseCrossover.");
            end
            if options.PopulationSize < 2
                error("radia:optuna:NSGAIIPopulation", ...
                    "PopulationSize must be at least 2.");
            end
            if options.PopulationSize < crossover.NParents
                error("radia:optuna:NSGAIIPopulation", ...
                    "PopulationSize must be at least the crossover parent count (%d).", ...
                    crossover.NParents);
            end
            obj.Seed = radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream = ...
                radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.IndependentSampler=radia.optuna.RandomSampler(options.Seed);
            obj.PopulationSize = options.PopulationSize;
            obj.MutationProbability = mutationProbability;
            obj.CrossoverProbability = options.CrossoverProbability;
            obj.SwappingProbability = options.SwappingProbability;
            obj.Crossover = crossover;
            obj.MutationScale = options.MutationScale;
            obj.ConstraintsFcn = options.ConstraintsFcn;
            obj.ElitePopulationSelectionStrategy = ...
                options.ElitePopulationSelectionStrategy;
            obj.ChildGenerationStrategy = options.ChildGenerationStrategy;
            obj.AfterTrialStrategy = options.AfterTrialStrategy;
        end

        function reseed_rng(obj)
            obj.IndependentSampler.reseed_rng();
            freshSeed=radia.optuna.internal.resolveSeed([]);
            obj.Stream=radia.optuna.internal.NumpyRandomState(freshSeed);
        end

        function searchSpace = inferRelativeSearchSpace(~,study,trial) %#ok<INUSD>
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.calculate( ...
                study,IncludePruned=false);
        end

        function searchSpace = infer_relative_search_space(obj,study,trial)
            if nargin < 3
                trial = [];
            end
            searchSpace = obj.inferRelativeSearchSpace(study,trial);
        end

        function beforeTrial(obj,study,trial)
            obj.IndependentSampler.beforeTrial(study,trial);
            obj.attach(study);
            generation = obj.assignGeneration(study,trial.Number);
            trial.setSystemAttr("nsgaii_generation",generation);
            trial.setSystemAttr("NSGAIISampler:generation",generation);
            if generation == 0
                obj.markFallback(trial,"initial_generation");
                obj.recordState(study,trial.Number);
                return
            end

            parentNumbers = obj.parentPoolForGeneration(study,generation);
            searchSpace = obj.inferRelativeSearchSpace(study,trial);
            if isempty(searchSpace) || ...
                    numel(parentNumbers) < obj.Crossover.NParents
                obj.markFallback(trial,"incompatible_parent_population");
                obj.recordState(study,trial.Number);
                return
            end
            [relativeSpace,child,mutatedNames,selectedParents] = ...
                obj.makeChild(study,parentNumbers,searchSpace);
            if ~isempty(relativeSpace)
                trial.setRelativeParameters(relativeSpace,child,"nsgaii");
            end
            trial.setSystemAttr("nsgaii_sampling_mode","joint");
            trial.setSystemAttr("nsgaii_joint_search_space", ...
                reshape([searchSpace.name],1,[]));
            trial.setSystemAttr("nsgaii_mutated_parameters",mutatedNames);
            trial.setSystemAttr("nsgaii_parent_trial_numbers", ...
                reshape(selectedParents,1,[]));
            obj.recordState(study,trial.Number);
        end

        function value = sampleFloat(obj,study,trial,name,low,high,options)
            value=obj.IndependentSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value = sampleInteger(obj,study,trial,name,low,high)
            value=obj.IndependentSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value = sampleCategorical(obj,study,trial,name,choices)
            value=obj.IndependentSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function values = sampleJoint(obj,study,trial,names,lows,highs,options)
            values = zeros(1,numel(names));
            for index = 1:numel(names)
                values(index)=obj.IndependentSampler.sampleFloat( ...
                    study,trial,names(index),lows(index),highs(index), ...
                    struct("Log",options.Log(index),"Step",NaN));
            end
        end

        function afterTrial(obj,study,trial)
            if ~isempty(obj.AfterTrialStrategy)
                obj.AfterTrialStrategy(study,trial,trial.State,trial.Values);
            elseif ismember(trial.State,["COMPLETE","PRUNED"]) && ...
                    ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial,obj.ConstraintsFcn(trial));
            end
            obj.recordState(study,trial.Number);
        end

        function population=select_parent(obj,study,generation)
            cache=obj.parentCacheIndex(generation);
            saved=[];
            if ~isempty(cache)
                saved=obj.ParentCaches(cache);
                obj.ParentCaches(cache)=[];
            end
            numbers=obj.parentPoolForGeneration(study,generation);
            generated=obj.parentCacheIndex(generation);
            if ~isempty(generated)
                obj.ParentCaches(generated)=[];
            end
            if ~isempty(saved)
                obj.ParentCaches(end+1)=saved;
            end
            trials=study.get_trials();
            indices=zeros(1,numel(numbers));
            for index=1:numel(numbers)
                indices(index)=find([trials.Number]==numbers(index),1);
            end
            population=trials(indices);
        end
    end

    methods (Access=private)
        function attach(obj,study)
            changed=isempty(obj.AttachedStudy) || ~isequal(obj.AttachedStudy,study);
            if changed
                obj.AttachedStudy=study;
                obj.Stream= ...
                    radia.optuna.internal.NumpyRandomState(obj.Seed);
                obj.GenerationTrialNumbers=zeros(0,1);
                obj.GenerationAssignments=zeros(0,1);
                obj.ParentCaches=struct("generation",{},"trial_numbers",{});
                obj.Restored=false;
            end
            if obj.Restored, return, end
            state=study.samplerState(obj.SamplerName,obj.StateSchema);
            if ~isempty(state)
                obj.restoreState(state);
            else
                state=study.samplerState(obj.SamplerName,obj.LegacyStateSchemaV2);
                if ~isempty(state)
                    obj.restoreLegacyStateV2(state);
                else
                    state=study.samplerState(obj.SamplerName,obj.LegacyStateSchemaV1);
                    if ~isempty(state), obj.restoreLegacyStateV1(state); end
                end
            end
            obj.reconcileExistingTrials(study);
            obj.Restored=true;
        end

        function restoreState(obj,state)
            required=["schema","seed","random_state","population_size", ...
                "mutation_probability","crossover_probability", ...
                "swapping_probability","crossover","strategies", ...
                "generation_by_trial","generation_parent_cache"];
            compatible=isstruct(state) && isscalar(state) && ...
                all(isfield(state,required));
            if compatible
                compatible=string(state.schema)==obj.StateSchema && ...
                    double(state.seed)==obj.Seed && ...
                    double(state.population_size)==obj.PopulationSize && ...
                    isequaln(double(state.mutation_probability),obj.MutationProbability) && ...
                    double(state.crossover_probability)==obj.CrossoverProbability && ...
                    double(state.swapping_probability)==obj.SwappingProbability && ...
                    isequaln(state.crossover,obj.Crossover.configuration()) && ...
                    isequaln(state.strategies,obj.strategyConfiguration());
            end
            if ~compatible
                error("radia:optuna:NSGAIIState", ...
                    "Stored NSGA-II sampler state is invalid or incompatible.");
            end
            [assignments,caches]=obj.validatePersistentState( ...
                state.generation_by_trial,state.generation_parent_cache);
            obj.Stream.State=state.random_state;
            obj.GenerationTrialNumbers=assignments(:,1);
            obj.GenerationAssignments=assignments(:,2);
            obj.ParentCaches=caches;
        end

        function restoreLegacyStateV2(obj,state)
            required=["schema","seed","random_state","population_size", ...
                "mutation_probability","crossover_probability", ...
                "generation_by_trial","generation_parent_cache"];
            compatible=isstruct(state) && isscalar(state) && ...
                all(isfield(state,required)) && ...
                string(state.schema)==obj.LegacyStateSchemaV2 && ...
                double(state.seed)==obj.Seed && ...
                double(state.population_size)==obj.PopulationSize && ...
                isequaln(double(state.mutation_probability),obj.MutationProbability) && ...
                double(state.crossover_probability)==obj.CrossoverProbability && ...
                isa(obj.Crossover,"radia.optuna.nsgaii.UniformCrossover") && ...
                obj.SwappingProbability==0.5 && obj.strategiesAreDefault();
            if ~compatible
                error("radia:optuna:NSGAIIState", ...
                    "Stored NSGA-II v2 state is incompatible.");
            end
            legacyCaches=state.generation_parent_cache;
            caches=struct("generation",{},"trial_numbers",{});
            for index=1:numel(legacyCaches)
                caches(end+1)=struct("generation",legacyCaches(index).generation, ...
                    "trial_numbers",legacyCaches(index).trial_numbers); %#ok<AGROW>
            end
            [assignments,caches]=obj.validatePersistentState( ...
                state.generation_by_trial,caches);
            obj.Stream.State=state.random_state;
            obj.GenerationTrialNumbers=assignments(:,1);
            obj.GenerationAssignments=assignments(:,2);
            obj.ParentCaches=caches;
        end

        function restoreLegacyStateV1(obj,state)
            if ~isstruct(state) || ~isscalar(state) || ...
                    ~all(isfield(state,["schema","seed","random_state"])) || ...
                    string(state.schema)~=obj.LegacyStateSchemaV1 || ...
                    double(state.seed)~=obj.Seed
                error("radia:optuna:NSGAIIState", ...
                    "Stored legacy NSGA-II state is invalid or incompatible.");
            end
            obj.Stream.State=state.random_state;
        end

        function [assignments,caches]=validatePersistentState(~,assignments,caches)
            assignments=double(assignments);
            if isempty(assignments), assignments=zeros(0,2); end
            if size(assignments,2)~=2 || any(~isfinite(assignments),"all") || ...
                    any(assignments~=floor(assignments),"all") || ...
                    any(assignments<0,"all") || ...
                    numel(unique(assignments(:,1)))~=size(assignments,1)
                error("radia:optuna:NSGAIIState", ...
                    "Stored generation assignments are invalid.");
            end
            if ~isstruct(caches) || (~isempty(caches) && ...
                    ~all(isfield(caches,["generation","trial_numbers"])))
                error("radia:optuna:NSGAIIState", ...
                    "Stored generation parent cache is invalid.");
            end
            normalized=struct("generation",{},"trial_numbers",{});
            for index=1:numel(caches)
                generation=double(caches(index).generation);
                numbers=reshape(double(caches(index).trial_numbers),[],1);
                if ~isscalar(generation) || ~isfinite(generation) || ...
                        generation<1 || generation~=floor(generation) || ...
                        any(~isfinite(numbers)) || any(numbers<0) || ...
                        any(numbers~=floor(numbers)) || ...
                        numel(unique(numbers))~=numel(numbers)
                    error("radia:optuna:NSGAIIState", ...
                        "Stored generation parent cache is invalid.");
                end
                normalized(end+1)=struct("generation",generation, ...
                    "trial_numbers",numbers); %#ok<AGROW>
            end
            caches=normalized;
            if ~isempty(caches) && ...
                    numel(unique([caches.generation]))~=numel(caches)
                error("radia:optuna:NSGAIIState", ...
                    "Stored parent-cache generations must be unique.");
            end
        end

        function config=strategyConfiguration(obj)
            config=struct( ...
                "constraints",obj.functionName(obj.ConstraintsFcn), ...
                "elite_population_selection", ...
                obj.functionName(obj.ElitePopulationSelectionStrategy), ...
                "child_generation",obj.functionName(obj.ChildGenerationStrategy), ...
                "after_trial",obj.functionName(obj.AfterTrialStrategy));
        end

        function result=strategiesAreDefault(obj)
            result=isempty(obj.ConstraintsFcn) && ...
                isempty(obj.ElitePopulationSelectionStrategy) && ...
                isempty(obj.ChildGenerationStrategy) && ...
                isempty(obj.AfterTrialStrategy);
        end

        function name=functionName(~,handle)
            if isempty(handle), name=""; else, name=string(func2str(handle)); end
        end

        function state=snapshot(obj)
            state=struct("schema",obj.StateSchema,"seed",obj.Seed, ...
                "population_size",obj.PopulationSize, ...
                "mutation_probability",obj.MutationProbability, ...
                "crossover_probability",obj.CrossoverProbability, ...
                "swapping_probability",obj.SwappingProbability, ...
                "crossover",obj.Crossover.configuration(), ...
                "strategies",obj.strategyConfiguration(), ...
                "random_state",obj.Stream.State, ...
                "generation_by_trial", ...
                [obj.GenerationTrialNumbers,obj.GenerationAssignments], ...
                "generation_parent_cache",obj.ParentCaches);
        end

        function recordState(obj,study,trialNumber)
            obj.attach(study);
            generation=0;
            index=find(obj.GenerationTrialNumbers==trialNumber,1);
            if ~isempty(index), generation=obj.GenerationAssignments(index); end
            study.recordSamplerState(obj.SamplerName,obj.StateSchema, ...
                trialNumber,generation,obj.snapshot());
        end

        function reconcileExistingTrials(obj,study)
            if isempty(study.TrialTable), return, end
            [trialNumbers,order]=sort(study.TrialTable.TrialNumber);
            states=study.TrialTable.State(order);
            if isempty(obj.GenerationTrialNumbers)
                generation=0;
                completed=0;
                for index=1:numel(trialNumbers)
                    obj.GenerationTrialNumbers(end+1,1)=trialNumbers(index);
                    obj.GenerationAssignments(end+1,1)=generation;
                    if states(index)=="COMPLETE"
                        completed=completed+1;
                        if completed>=obj.PopulationSize
                            generation=generation+1;
                            completed=0;
                        end
                    end
                end
                return
            end
            for number=reshape(trialNumbers,1,[])
                if ~ismember(number,obj.GenerationTrialNumbers)
                    obj.GenerationTrialNumbers(end+1,1)=number;
                    obj.GenerationAssignments(end+1,1)=obj.currentGeneration(study);
                end
            end
        end

        function generation=assignGeneration(obj,study,trialNumber)
            index=find(obj.GenerationTrialNumbers==trialNumber,1);
            if ~isempty(index)
                generation=obj.GenerationAssignments(index);
                return
            end
            generation=obj.currentGeneration(study);
            obj.GenerationTrialNumbers(end+1,1)=trialNumber;
            obj.GenerationAssignments(end+1,1)=generation;
        end

        function generation=currentGeneration(obj,study)
            if isempty(obj.GenerationAssignments), generation=0; return, end
            generation=max(obj.GenerationAssignments);
            numbers=obj.GenerationTrialNumbers( ...
                obj.GenerationAssignments==generation);
            complete=study.TrialTable.TrialNumber( ...
                study.TrialTable.State=="COMPLETE");
            if sum(ismember(numbers,complete))>=obj.PopulationSize
                generation=generation+1;
            end
        end

        function parentNumbers=parentPoolForGeneration(obj,study,generation)
            cache=obj.parentCacheIndex(generation);
            if ~isempty(cache)
                cached=reshape(obj.ParentCaches(cache).trial_numbers,[],1);
                parentNumbers=study.TrialTable.TrialNumber(ismember( ...
                    study.TrialTable.TrialNumber,cached));
                return
            end
            previous=generation-1;
            complete=study.TrialTable.TrialNumber( ...
                study.TrialTable.State=="COMPLETE");
            offspring=obj.GenerationTrialNumbers( ...
                obj.GenerationAssignments==previous & ...
                ismember(obj.GenerationTrialNumbers,complete));
            inherited=zeros(0,1);
            priorCache=obj.parentCacheIndex(previous);
            if ~isempty(priorCache)
                cached=reshape( ...
                    obj.ParentCaches(priorCache).trial_numbers,[],1);
                inherited=study.TrialTable.TrialNumber(ismember( ...
                    study.TrialTable.TrialNumber,cached));
            end
            % BaseGASampler.select_parent passes the just-finished
            % generation before the inherited parent population.  The
            % order is observable when ranks/crowding ties are stable and
            % later tournament draws index that parent list.
            candidates=unique([offspring;inherited],"stable");
            values=obj.objectivesForTrials(study,candidates);
            valid=all(isfinite(values),2);
            candidates=candidates(valid);
            values=values(valid,:);
            if isempty(obj.ElitePopulationSelectionStrategy)
                order=radia.optuna.internal.ParetoSupport.eliteSelectionOrder( ...
                    study,candidates,values,obj.PopulationSize);
                parentNumbers=candidates(order);
            else
                parentNumbers=reshape(double( ...
                    obj.ElitePopulationSelectionStrategy( ...
                    study,candidates)),[],1);
                if numel(parentNumbers)>obj.PopulationSize || ...
                        numel(unique(parentNumbers))~=numel(parentNumbers) || ...
                        any(~ismember(parentNumbers,candidates))
                    error("radia:optuna:NSGAIIStrategy", ...
                        "Elite selection returned an invalid parent population.");
                end
            end
            obj.ParentCaches(end+1)=struct("generation",generation, ...
                "trial_numbers",reshape(parentNumbers,[],1));
        end

        function index=parentCacheIndex(obj,generation)
            if isempty(obj.ParentCaches), index=[]; return, end
            index=find([obj.ParentCaches.generation]==generation,1);
        end

        function values=objectivesForTrials(~,study,trialNumbers)
            values=NaN(numel(trialNumbers),numel(study.Directions));
            for row=1:numel(trialNumbers)
                for objective=1:numel(study.Directions)
                    selected=study.ObjectiveTable.TrialNumber==trialNumbers(row) & ...
                        study.ObjectiveTable.ObjectiveIndex==objective;
                    if sum(selected)==1
                        values(row,objective)=study.ObjectiveTable.Value(selected);
                    end
                end
            end
        end

        function [relativeSpace,child,mutatedNames,selectedParents]= ...
                makeChild(obj,study,parentNumbers,searchSpace)
            if ~isempty(obj.ChildGenerationStrategy)
                child=obj.ChildGenerationStrategy( ...
                    study,searchSpace,parentNumbers);
                if isstruct(child) && isscalar(child)
                    converted=cell(1,numel(searchSpace));
                    for index=1:numel(searchSpace)
                        key=matlab.lang.makeValidName(searchSpace(index).name);
                        if ~isfield(child,key)
                            error("radia:optuna:NSGAIIStrategy", ...
                                "Child strategy omitted parameter '%s'.", ...
                                searchSpace(index).name);
                        end
                        converted{index}=child.(key);
                    end
                    child=converted;
                end
                relativeSpace=searchSpace;
                mutatedNames=strings(1,0);
                selectedParents=reshape(parentNumbers,1,[]);
                obj.validateChild(relativeSpace,child,searchSpace);
                return
            end
            if rand(obj.Stream)<obj.CrossoverProbability
                usedCrossover=true;
                [child,selectedParents]=obj.performCrossover( ...
                    study,parentNumbers,searchSpace);
            else
                usedCrossover=false;
                selectedParents=parentNumbers( ...
                    randi(obj.Stream,numel(parentNumbers)));
                child=obj.parameterValues(study,selectedParents,searchSpace);
            end
            dimension=numel(searchSpace);
            probability=obj.MutationProbability;
            if isnan(probability), probability=1/max(1,dimension); end
            % Optuna constructs the crossover result by inserting the
            % categorical parameters first and the numerical parameters
            % second.  Python preserves that dict insertion order when it
            % consumes one mutation draw per parameter.  Keep the public
            % search-space order, but assign the draws in the same order so
            % a fixed seed follows the upstream proposal sequence exactly.
            if usedCrossover
                categorical=arrayfun(@(x) ...
                    x.distribution.kind=="categorical",searchSpace);
                mutationOrder=[reshape(find(categorical),1,[]), ...
                    reshape(find(~categorical),1,[])];
            else
                % The no-crossover branch copies a dict comprehension in
                % the original search-space order.
                mutationOrder=1:dimension;
            end
            mutationDraws=rand(obj.Stream,1,dimension)<probability;
            mutated=false(1,dimension);
            mutated(mutationOrder)=mutationDraws;
            mutatedNames=reshape([searchSpace(mutated).name],1,[]);
            relativeSpace=searchSpace(~mutated);
            child=child(~mutated);
        end

        function [child,selectedParents]=performCrossover( ...
                obj,study,parentNumbers,searchSpace)
            for attempt=1:10000
                selectedParents=obj.selectParents( ...
                    study,parentNumbers,obj.Crossover.NParents);
                parents=obj.parameterValues(study,selectedParents,searchSpace);
                child=cell(1,numel(searchSpace));
                categorical=arrayfun(@(x) ...
                    x.distribution.kind=="categorical",searchSpace);
                if any(categorical)
                    useLast=rand(obj.Stream,1,sum(categorical))>= ...
                        obj.SwappingProbability;
                    positions=find(categorical);
                    for offset=1:numel(positions)
                        row=1+double(useLast(offset));
                        child{positions(offset)}=parents{row,positions(offset)};
                    end
                end
                numeric=find(~categorical);
                if ~isempty(numeric)
                    [transformed,bounds]=obj.transformParents( ...
                        parents(:,numeric),searchSpace(numeric));
                    crossed=obj.Crossover.crossover( ...
                        transformed,obj.Stream,study,bounds);
                    if numel(crossed)~=numel(numeric) || ...
                            any(~isfinite(crossed))
                        continue
                    end
                    numericChild=obj.untransformChild( ...
                        reshape(crossed,1,[]),searchSpace(numeric));
                    if isempty(numericChild), continue, end
                    child(numeric)=numericChild;
                end
                if obj.childContained(child,searchSpace), return, end
            end
            error("radia:optuna:NSGAIICrossover", ...
                "Unable to generate a child inside the search space.");
        end

        function selected=selectParents(obj,study,population,count)
            selected=zeros(0,1);
            values=obj.objectivesForTrials(study,population);
            for parent=1:count
                available=population(~ismember(population,selected));
                candidateIndices=randi(obj.Stream,numel(available),1,2);
                left=available(candidateIndices(1));
                right=available(candidateIndices(2));
                leftRow=find(population==left,1);
                rightRow=find(population==right,1);
                if radia.optuna.internal.ParetoSupport.constrainedDominates( ...
                        study,left,values(leftRow,:),right,values(rightRow,:))
                    selected(end+1,1)=left; %#ok<AGROW>
                else
                    selected(end+1,1)=right; %#ok<AGROW>
                end
            end
        end

        function population=parameterValues(obj,study,trialNumbers,searchSpace)
            population=cell(numel(trialNumbers),numel(searchSpace));
            for row=1:numel(trialNumbers)
                for column=1:numel(searchSpace)
                    selected=study.ParamTable.TrialNumber==trialNumbers(row) & ...
                        study.ParamTable.Name==searchSpace(column).name;
                    if sum(selected)~=1
                        error("radia:optuna:NSGAIIObservations", ...
                            "Parent trial %d lacks parameter '%s'.", ...
                            trialNumbers(row),searchSpace(column).name);
                    end
                    tableRow=study.ParamTable(selected,:);
                    distribution=radia.optuna.internal.DistributionCodec.decode( ...
                        tableRow.Kind,tableRow.Distribution);
                    if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                            distribution,searchSpace(column).distribution)
                        error("radia:optuna:NSGAIIObservations", ...
                            "Parent trial %d has an incompatible distribution.", ...
                            trialNumbers(row));
                    end
                    [accepted,value]=obj.parameterValue(tableRow,distribution);
                    if ~accepted
                        error("radia:optuna:NSGAIIObservations", ...
                            "Parent trial %d contains an invalid parameter.", ...
                            trialNumbers(row));
                    end
                    population{row,column}=value;
                end
            end
        end

        function [accepted,value]=parameterValue(~,row,distribution)
            accepted=false;
            if distribution.kind=="categorical"
                if isfinite(row.ValueNumeric)
                    token="numeric:"+string(row.ValueNumeric,"%.17g");
                else
                    token=row.ValueText;
                end
                tokens=radia.optuna.internal.DistributionCodec.choiceTokens( ...
                    distribution.choices);
                index=find(tokens==token,1);
                if isempty(index), value=[]; return, end
                value=radia.optuna.internal.DistributionCodec.choiceAt( ...
                    distribution.choices,index);
                accepted=true;
                return
            end
            value=row.ValueNumeric;
            accepted=isfinite(value) && value>=distribution.low && ...
                value<=distribution.high && (~distribution.log || value>0);
            if accepted && isfinite(distribution.step)
                grid=(value-distribution.low)/distribution.step;
                accepted=abs(grid-round(grid))<=1e-10*max(1,abs(grid));
            end
            if accepted && distribution.kind=="integer"
                accepted=value==round(value);
                value=round(value);
            end
        end

        function [parents,bounds]=transformParents(~,population,searchSpace)
            parents=zeros(size(population));
            bounds=zeros(numel(searchSpace),2);
            for column=1:numel(searchSpace)
                distribution=searchSpace(column).distribution;
                values=cellfun(@double,population(:,column));
                if distribution.log, values=log(values); end
                parents(:,column)=values;
                low=distribution.low;
                high=distribution.high;
                step=distribution.step;
                if distribution.kind=="integer" && distribution.log
                    half=0.5*step;
                    bounds(column,:)=[log(low-half),log(high+half)];
                else
                    if distribution.log
                        low=log(low); high=log(high);
                    end
                    half=0;
                    if isfinite(step), half=0.5*step; end
                    bounds(column,:)=[low-half,high+half];
                end
            end
        end

        function child=untransformChild(obj,transformed,searchSpace)
            child=cell(1,numel(searchSpace));
            for index=1:numel(searchSpace)
                distribution=searchSpace(index).distribution;
                value=transformed(index);
                if value<obj.transformedBound(distribution,1) || ...
                        value>obj.transformedBound(distribution,2)
                    child={}; return
                end
                if distribution.log, value=exp(value); end
                if distribution.kind=="float"
                    if isfinite(distribution.step)
                        value=distribution.low+obj.roundTiesToEven( ...
                            (value-distribution.low)/distribution.step)* ...
                            distribution.step;
                        value=min(max(value,distribution.low),distribution.high);
                    elseif distribution.low~=distribution.high
                        value=min(value,obj.nextDown(distribution.high));
                    end
                elseif distribution.log
                    value=obj.roundTiesToEven(value);
                    value=min(max(value,distribution.low),distribution.high);
                else
                    value=distribution.low+obj.roundTiesToEven( ...
                        (value-distribution.low)/distribution.step)* ...
                        distribution.step;
                    value=min(max(value,distribution.low),distribution.high);
                    value=round(value);
                end
                child{index}=value;
            end
        end

        function bound=transformedBound(~,distribution,side)
            half=0;
            if isfinite(distribution.step), half=0.5*distribution.step; end
            if side==1
                value=distribution.low;
                offset=-half;
            else
                value=distribution.high;
                offset=half;
            end
            if distribution.kind=="integer" && distribution.log
                bound=log(value+offset);
            elseif distribution.log
                bound=log(value)+offset;
            else
                bound=value+offset;
            end
        end

        function result=childContained(~,child,searchSpace)
            result=numel(child)==numel(searchSpace);
            for index=1:numel(searchSpace)
                if ~result, return, end
                distribution=searchSpace(index).distribution;
                value=child{index};
                if distribution.kind=="categorical"
                    token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                    result=ismember(token, ...
                        radia.optuna.internal.DistributionCodec.choiceTokens( ...
                        distribution.choices));
                else
                    result=isnumeric(value) && isscalar(value) && isfinite(value) && ...
                        value>=distribution.low && value<=distribution.high;
                end
            end
        end

        function validateChild(~,relativeSpace,child,searchSpace)
            if ~isstruct(relativeSpace) || ~iscell(child) || ...
                    numel(relativeSpace)~=numel(child) || ...
                    any(~ismember([relativeSpace.name],[searchSpace.name]))
                error("radia:optuna:NSGAIIStrategy", ...
                    "Child-generation strategy returned an invalid child.");
            end
        end

        function value=randomValue(obj,distribution)
            if distribution.kind=="categorical"
                index=randi(obj.Stream,numel(distribution.choices));
                value=radia.optuna.internal.DistributionCodec.choiceAt( ...
                    distribution.choices,index);
                return
            end
            if distribution.low==distribution.high
                value=distribution.low;
            else
                value=obj.uniform(distribution.low,distribution.high, ...
                    distribution.log);
            end
            value=obj.quantize(value,distribution.low, ...
                distribution.high,distribution.step);
            if distribution.kind=="integer", value=round(value); end
        end

        function value=uniform(obj,low,high,logScale)
            fraction=rand(obj.Stream);
            if logScale
                value=exp(log(low)+fraction*(log(high)-log(low)));
            else
                value=low+fraction*(high-low);
            end
        end

        function value=quantize(obj,value,low,high,step)
            if isfinite(step)
                value=low+obj.roundTiesToEven((value-low)/step)*step;
            end
            value=min(max(value,low),high);
        end

        function value=roundTiesToEven(~,value)
            value=radia.optuna.internal.UpstreamNumerics.roundTiesToEven( ...
                value);
        end

        function value=nextDown(~,value)
            bits=typecast(double(value),'uint64');
            if value>0
                bits=bits-uint64(1);
            elseif value<0
                bits=bits+uint64(1);
            else
                bits=bitor(bitshift(uint64(1),63),uint64(1));
            end
            value=typecast(bits,'double');
        end

        function markFallback(~,trial,reason)
            trial.setSystemAttr("nsgaii_sampling_mode","independent_fallback");
            trial.setSystemAttr("nsgaii_fallback_reason",reason);
        end
    end
end
