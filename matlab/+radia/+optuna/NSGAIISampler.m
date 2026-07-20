classdef NSGAIISampler < handle
    %NSGAIISAMPLER Pareto-ranked elitist crossover and mutation sampler.
    properties (SetAccess=private)
        Stream
        PopulationSize (1,1) double = 24
        MutationProbability (1,1) double = 0.15
        CrossoverProbability (1,1) double = 0.9
        MutationScale (1,1) double = 0.1
    end
    methods
        function obj=NSGAIISampler(options)
            arguments
                options.Seed (1,1) double=0
                options.PopulationSize (1,1) double {mustBeInteger,mustBePositive}=24
                options.MutationProbability (1,1) double=0.15
                options.CrossoverProbability (1,1) double=0.9
                options.MutationScale (1,1) double {mustBePositive}=0.1
            end
            if options.MutationProbability<0 || options.MutationProbability>1 || ...
                    options.CrossoverProbability<0 || options.CrossoverProbability>1
                error("radia:optuna:NSGAIIProbability","Mutation and crossover probabilities must be in [0,1].");
            end
            obj.Stream=RandStream("mt19937ar","Seed",options.Seed); obj.PopulationSize=options.PopulationSize;
            obj.MutationProbability=options.MutationProbability; obj.CrossoverProbability=options.CrossoverProbability;
            obj.MutationScale=options.MutationScale;
        end
        function value=sampleFloat(obj,study,trial,name,low,high,options) %#ok<INUSD>
            obj.validate(low,high,options); [x,objectives]=radia.optuna.internal.ParetoSupport.numericObservations(study,name);
            if numel(x)<obj.PopulationSize
                value=obj.uniform(low,high,options.Log); value=obj.quantize(value,low,high,options.Step); return
            end
            order=radia.optuna.internal.ParetoSupport.preferenceOrder(objectives,study.Directions);
            elite=order(1:min(obj.PopulationSize,numel(order))); x=x(elite); objectives=objectives(elite,:);
            [rank,crowding]=radia.optuna.internal.ParetoSupport.rankAndCrowding(objectives,study.Directions);
            first=obj.tournament(rank,crowding); second=obj.tournament(rank,crowding);
            transform=@(v)v; inverse=@(v)v; lo=low; hi=high;
            if options.Log, transform=@log; inverse=@exp; lo=log(low); hi=log(high); end
            parents=transform([x(first),x(second)]);
            if rand(obj.Stream)<obj.CrossoverProbability
                alpha=rand(obj.Stream); child=alpha*parents(1)+(1-alpha)*parents(2);
            else
                child=parents(1);
            end
            if rand(obj.Stream)<obj.MutationProbability
                child=child+obj.MutationScale*(hi-lo)*randn(obj.Stream);
            end
            value=obj.quantize(inverse(min(max(child,lo),hi)),low,high,options.Step);
        end
        function value=sampleInteger(obj,study,trial,name,low,high)
            value=round(obj.sampleFloat(study,trial,name,low,high,struct("Log",false,"Step",1)));
        end
        function value=sampleCategorical(obj,study,trial,name,choices) %#ok<INUSD>
            if isempty(choices), error("radia:optuna:Choices","Choices must not be empty."); end
            [tokens,objectives]=radia.optuna.internal.ParetoSupport.categoricalObservations(study,name);
            if numel(tokens)<obj.PopulationSize || rand(obj.Stream)<obj.MutationProbability
                value=obj.choiceAt(choices,1+floor(rand(obj.Stream)*numel(choices))); return
            end
            order=radia.optuna.internal.ParetoSupport.preferenceOrder(objectives,study.Directions);
            elite=order(1:min(obj.PopulationSize,numel(order))); tokens=tokens(elite); objectives=objectives(elite,:);
            [rank,crowding]=radia.optuna.internal.ParetoSupport.rankAndCrowding(objectives,study.Directions);
            token=tokens(obj.tournament(rank,crowding));
            for k=1:numel(choices)
                if obj.token(obj.choiceAt(choices,k))==token, value=obj.choiceAt(choices,k); return, end
            end
            value=obj.choiceAt(choices,1+floor(rand(obj.Stream)*numel(choices)));
        end
        function beforeTrial(obj,study,trial) %#ok<INUSD>
        end
        function afterTrial(obj,study,trial) %#ok<INUSD>
        end
    end
    methods (Access=private)
        function index=tournament(obj,rank,crowding)
            candidates=1+floor(rand(obj.Stream,1,2)*numel(rank)); a=candidates(1); b=candidates(2);
            if rank(a)<rank(b) || (rank(a)==rank(b) && crowding(a)>=crowding(b)), index=a; else, index=b; end
        end
        function validate(~,low,high,options)
            if ~(isfinite(low)&&isfinite(high)&&low<high), error("radia:optuna:Bounds","Float bounds must satisfy low < high."); end
            if options.Log&&low<=0, error("radia:optuna:LogBounds","Log bounds must be positive."); end
            if isfinite(options.Step)&&options.Step<=0, error("radia:optuna:Step","Step must be positive."); end
        end
        function value=uniform(obj,low,high,logScale)
            u=rand(obj.Stream); if logScale, value=exp(log(low)+u*(log(high)-log(low))); else, value=low+u*(high-low); end
        end
        function value=quantize(~,value,low,high,step)
            if isfinite(step), value=low+round((value-low)/step)*step; end, value=min(max(value,low),high);
        end
        function token=token(~,value)
            if isnumeric(value)&&isscalar(value), token="numeric:"+string(value,"%.17g"); else, token=string(jsonencode(value)); end
        end
        function value=choiceAt(~,choices,index)
            if iscell(choices), value=choices{index}; else, value=choices(index); end
        end
    end
end
