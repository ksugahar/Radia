classdef NSGAIIISampler < handle
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
                options.Seed (1,1) double = 0
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
            obj.Seed=double(options.Seed);
            obj.PopulationSize=double(options.PopulationSize);
            obj.DividingParameter=double(options.DividingParameter);
            obj.ReferencePoints=referencePoints;
            obj.Stream=RandStream("mt19937ar","Seed",obj.Seed);
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
                if numel(front)<=remaining
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
            while numel(chosen)<count && ~isempty(available)
                active=unique(associations(available));
                minimum=min(nicheCounts(active));
                tied=active(nicheCounts(active)==minimum);
                reference=tied(randi(obj.Stream,numel(tied)));
                local=available(associations(available)==reference);
                if nicheCounts(reference)==0
                    direction=references(reference,:);
                    projection=(points(local,:)*direction')*direction;
                    distance=vecnorm(points(local,:)-projection,2,2);
                    [~,pick]=min(distance);
                else
                    pick=randi(obj.Stream,numel(local));
                end
                selectedIndex=local(pick);
                chosen(end+1,1)=selectedIndex; %#ok<AGROW>
                available(available==selectedIndex)=[];
                nicheCounts(reference)=nicheCounts(reference)+1;
            end
        end
    end
end
