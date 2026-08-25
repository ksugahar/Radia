classdef NSGAIIISampler < radia.optuna.BaseSampler
    %NSGAIIISAMPLER Constrained joint NSGA-III with reference-line niching.
    %   Child generation, dynamic-space fallback, generation caching, and
    %   constrained Pareto ranks are shared with NSGAIISampler. The elite
    %   cutoff front is selected by Optuna-style NSGA-III normalization and
    %   reference-line niche preservation.

    properties (SetAccess=private)
        Seed (1,1) double = 0
        PopulationSize (1,1) double = 50
        DividingParameter (1,1) double = 3
        ReferencePoints double = zeros(0,0)
        Core (1,1) radia.optuna.NSGAIISampler
    end

    properties (Access=private)
        Stream
    end

    methods
        function obj=NSGAIIISampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.PopulationSize (1,1) double ...
                    {mustBeInteger,mustBePositive} = 50
                options.DividingParameter (1,1) double ...
                    {mustBeInteger,mustBePositive} = 3
                options.ReferencePoints double = zeros(0,0)
                options.MutationProbability (1,1) double = NaN
                options.Crossover = []
                options.CrossoverProbability (1,1) double = 0.9
                options.SwappingProbability (1,1) double = 0.5
                options.ConstraintsFcn = []
                options.ChildGenerationStrategy = []
                options.AfterTrialStrategy = []
            end
            if options.PopulationSize<2
                error("radia:optuna:NSGAIIIPopulation", ...
                    "PopulationSize must be at least 2.");
            end
            referencePoints=double(options.ReferencePoints);
            if ~isempty(referencePoints) && ...
                    (ismatrix(referencePoints)==false || ...
                    any(~isfinite(referencePoints),"all") || ...
                    any(referencePoints<0,"all") || ...
                    any(sum(referencePoints,2)<=0))
                error("radia:optuna:NSGAIIIReferencePoints", ...
                    "ReferencePoints must be finite, nonnegative, nonzero rows.");
            end
            obj.Seed=radia.optuna.internal.resolveSeed(options.Seed);
            obj.PopulationSize=double(options.PopulationSize);
            obj.DividingParameter=double(options.DividingParameter);
            obj.ReferencePoints=referencePoints;
            obj.Stream=radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.Core=radia.optuna.NSGAIISampler( ...
                Seed=obj.Seed,PopulationSize=obj.PopulationSize, ...
                MutationProbability=options.MutationProbability, ...
                Crossover=options.Crossover, ...
                CrossoverProbability=options.CrossoverProbability, ...
                SwappingProbability=options.SwappingProbability, ...
                ConstraintsFcn=options.ConstraintsFcn, ...
                ElitePopulationSelectionStrategy=@obj.selectElite, ...
                ChildGenerationStrategy=options.ChildGenerationStrategy, ...
                AfterTrialStrategy=options.AfterTrialStrategy);
            % Use one stream for niche tie-breaking and child generation so
            % the persisted NSGA-II core state also restores NSGA-III.
            obj.Stream=obj.Core.Stream;
        end

        function searchSpace=inferRelativeSearchSpace(obj,study,trial)
            searchSpace=obj.Core.inferRelativeSearchSpace(study,trial);
        end

        function searchSpace=infer_relative_search_space(obj,study,trial)
            if nargin<3, trial=[]; end
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
        end

        function beforeTrial(obj,study,trial)
            obj.Core.beforeTrial(study,trial);
            attrs=trial.SystemAttrs;
            if isfield(attrs,"nsgaii_generation")
                trial.setSystemAttr("nsgaiii_generation", ...
                    attrs.nsgaii_generation);
            end
            trial.setSystemAttr("nsgaiii_elite_strategy", ...
                "reference_line_niching");
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options)
            value=obj.Core.sampleFloat(study,trial,name,low,high,options);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            value=obj.Core.sampleInteger(study,trial,name,low,high);
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            value=obj.Core.sampleCategorical(study,trial,name,choices);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            values=obj.Core.sampleJoint(study,trial,names,lows,highs,options);
        end

        function afterTrial(obj,study,trial)
            obj.Core.afterTrial(study,trial);
        end

        function selectedNumbers=selectElitePopulation(obj,study,candidateNumbers)
            %SELECTELITEPOPULATION Expose the deterministic elite contract.
            selectedNumbers=obj.selectElite(study,candidateNumbers);
        end
    end

    methods (Access=private)
        function selectedNumbers=selectElite(obj,study,candidateNumbers)
            % Core.attach may replace its RandStream while restoring a study.
            % Rebind here before any niche tie-break so replay is persisted.
            obj.Stream=obj.Core.Stream;
            candidateNumbers=reshape(double(candidateNumbers),[],1);
            values=obj.objectives(study,candidateNumbers);
            [~,~,rank]=radia.optuna.internal.ParetoSupport. ...
                constrainedRankAndCrowding(study,candidateNumbers,values);
            selected=zeros(0,1);
            cutoff=zeros(0,1);
            for level=reshape(unique(rank,"sorted"),1,[])
                front=find(rank==level);
                remaining=obj.PopulationSize-numel(selected);
                if remaining<=0, break, end
                % Upstream runs reference-line preservation even when the
                % cutoff front exactly fills the remaining population.  In
                % addition to ordering the parents, that operation consumes
                % the seeded RNG through its bucket shuffles.
                if numel(front)<remaining
                    selected=[selected;front]; %#ok<AGROW>
                else
                    cutoff=front;
                    break
                end
            end
            if ~isempty(cutoff)
                referencePoints=obj.referencePointsFor( ...
                    numel(study.Directions));
                % Optuna normalizes only the fronts that can still affect
                % survival: the already accepted fronts plus the cutoff
                % front.  Worse fronts must not move the ideal point or
                % intercept plane.
                relevant=[selected;cutoff];
                normalized=radia.optuna.internal.NSGAIIISupport. ...
                    normalizeObjectives(values(relevant,:));
                associations=radia.optuna.internal.NSGAIIISupport. ...
                    associate(normalized,referencePoints);
                accepted=(1:numel(selected))';
                available=(numel(selected)+(1:numel(cutoff)))';
                chosen=obj.nicheSelect(associations,normalized, ...
                    referencePoints,accepted,available, ...
                    obj.PopulationSize-numel(selected));
                selected=[selected;relevant(chosen)];
            end
            selectedNumbers=candidateNumbers(selected);
        end

        function values=objectives(~,study,trialNumbers)
            values=NaN(numel(trialNumbers),numel(study.Directions));
            for row=1:numel(trialNumbers)
                for objective=1:numel(study.Directions)
                    mask=study.ObjectiveTable.TrialNumber==trialNumbers(row) & ...
                        study.ObjectiveTable.ObjectiveIndex==objective;
                    if sum(mask)==1
                        values(row,objective)=study.ObjectiveTable.Value(mask);
                    end
                end
            end
            if any(isnan(values),"all")
                error("radia:optuna:NSGAIIIObservations", ...
                    "NSGA-III objective values must not be NaN.");
            end
        end

        function points=referencePointsFor(obj,nObjectives)
            if ~isempty(obj.ReferencePoints)
                if size(obj.ReferencePoints,2)~=nObjectives
                    error("radia:optuna:NSGAIIIReferencePoints", ...
                        "ReferencePoints must have one column per objective.");
                end
                points=obj.ReferencePoints;
                return
            end
            points=radia.optuna.internal.NSGAIIISupport. ...
                defaultReferencePoints(nObjectives,obj.DividingParameter);
        end

        function chosen=nicheSelect(obj,associations,points,references, ...
                selected,available,count)
            references=references./vecnorm(references,2,2);
            nicheCounts=accumarray(associations(selected),1, ...
                [size(references,1),1]);
            chosen=zeros(0,1);
            referenceCount=size(references,1);
            borderline=cell(referenceCount,1);
            borderlineShuffled=false(referenceCount,1);
            for localIndex=1:numel(available)
                populationIndex=available(localIndex);
                reference=associations(populationIndex);
                direction=references(reference,:);
                projection=(points(populationIndex,:)*direction')*direction;
                distance=norm(points(populationIndex,:)-projection);
                borderline{reference}(end+1,:)=[distance,localIndex];
            end

            % Python defaultdict/list insertion order is the first
            % appearance of a reference in the cutoff population.
            activeReferences=unique(associations(available),"stable");
            buckets=cell(max(count+1,2),1);
            for reference=reshape(activeReferences,1,[])
                bucket=nicheCounts(reference)+1;
                if bucket>numel(buckets), buckets{bucket}=zeros(1,0); end
                buckets{bucket}(end+1)=reference;
            end

            bucketCount=-1;
            while numel(chosen)<count
                bucketIndex=bucketCount+1;
                if bucketIndex<1 || bucketIndex>numel(buckets) || ...
                        isempty(buckets{bucketIndex})
                    bucketCount=bucketCount+1;
                    bucketIndex=bucketCount+1;
                    if bucketIndex>numel(buckets)
                        buckets{bucketIndex}=zeros(1,0);
                    end
                    order=randperm(obj.Stream,numel(buckets{bucketIndex}));
                    buckets{bucketIndex}=buckets{bucketIndex}(order);
                    continue
                end

                % list.pop() removes the last shuffled reference.
                reference=buckets{bucketIndex}(end);
                buckets{bucketIndex}(end)=[];
                candidates=borderline{reference};
                if bucketCount==0
                    % sort(reverse=True); pop() chooses the smallest
                    % (distance, original-population-index) tuple.
                    [~,order]=sortrows(candidates,[1 2],"ascend");
                    pick=order(1);
                else
                    if ~borderlineShuffled(reference)
                        order=randperm(obj.Stream,size(candidates,1));
                        candidates=candidates(order,:);
                        borderline{reference}=candidates;
                        borderlineShuffled(reference)=true;
                    end
                    pick=size(candidates,1);
                end
                localIndex=candidates(pick,2);
                chosen(end+1,1)=available(localIndex); %#ok<AGROW>
                candidates(pick,:)=[];
                borderline{reference}=candidates;
                if ~isempty(candidates)
                    nextBucket=bucketCount+2;
                    if nextBucket>numel(buckets)
                        buckets{nextBucket}=zeros(1,0);
                    end
                    buckets{nextBucket}(end+1)=reference;
                end
            end
        end
    end
end
