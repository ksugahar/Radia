classdef MOTPESampler < handle
    %MOTPESAMPLER Multi-objective TPE using Pareto rank and crowding.
    properties (SetAccess=private)
        Stream
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.25
        NumberOfEIChoices (1,1) double = 32
        PriorWeight (1,1) double = 1
    end
    methods
        function obj=MOTPESampler(options)
            arguments
                options.Seed (1,1) double=0
                options.NStartupTrials (1,1) double {mustBeInteger,mustBeNonnegative}=10
                options.Gamma (1,1) double=0.25
                options.NumberOfEIChoices (1,1) double {mustBeInteger,mustBePositive}=32
                options.PriorWeight (1,1) double {mustBePositive}=1
            end
            if options.Gamma<=0 || options.Gamma>=1, error("radia:optuna:MOTPEGamma","Gamma must be between zero and one."); end
            obj.Stream=RandStream("mt19937ar","Seed",options.Seed); obj.NStartupTrials=options.NStartupTrials;
            obj.Gamma=options.Gamma; obj.NumberOfEIChoices=options.NumberOfEIChoices; obj.PriorWeight=options.PriorWeight;
        end
        function value=sampleFloat(obj,study,trial,name,low,high,options) %#ok<INUSD>
            obj.validate(low,high,options); [x,objectives]=radia.optuna.internal.ParetoSupport.numericObservations(study,name);
            if numel(x)<max(2,obj.NStartupTrials), value=obj.uniform(low,high,options.Log); value=obj.quantize(value,low,high,options.Step); return, end
            order=radia.optuna.internal.ParetoSupport.preferenceOrder(objectives,study.Directions);
            nGood=max(1,min(numel(x)-1,ceil(obj.Gamma*numel(x)))); good=x(order(1:nGood)); bad=x(order(nGood+1:end));
            transform=@(v)v; inverse=@(v)v; lo=low; hi=high;
            if options.Log, transform=@log; inverse=@exp; lo=log(low); hi=log(high); end
            good=transform(good); bad=transform(bad); bwGood=obj.bandwidth(good,hi-lo); bwBad=obj.bandwidth(bad,hi-lo);
            candidates=[good(:);lo+(hi-lo)*rand(obj.Stream,obj.NumberOfEIChoices,1)]; score=zeros(size(candidates));
            for k=1:numel(candidates), score(k)=obj.logKde(candidates(k),good,bwGood)-obj.logKde(candidates(k),bad,bwBad); end
            [~,best]=max(score); value=obj.quantize(inverse(candidates(best)),low,high,options.Step);
        end
        function value=sampleInteger(obj,study,trial,name,low,high)
            value=round(obj.sampleFloat(study,trial,name,low,high,struct("Log",false,"Step",1)));
        end
        function value=sampleCategorical(obj,study,trial,name,choices) %#ok<INUSD>
            if isempty(choices), error("radia:optuna:Choices","Choices must not be empty."); end
            [tokens,objectives]=radia.optuna.internal.ParetoSupport.categoricalObservations(study,name); count=numel(choices);
            if numel(tokens)<max(2,obj.NStartupTrials), value=obj.choiceAt(choices,1+floor(rand(obj.Stream)*count)); return, end
            order=radia.optuna.internal.ParetoSupport.preferenceOrder(objectives,study.Directions);
            nGood=max(1,min(numel(tokens)-1,ceil(obj.Gamma*numel(tokens)))); good=tokens(order(1:nGood)); bad=tokens(order(nGood+1:end)); scores=zeros(count,1);
            for k=1:count
                token=obj.token(obj.choiceAt(choices,k)); pg=(sum(good==token)+obj.PriorWeight)/(numel(good)+obj.PriorWeight*count);
                pb=(sum(bad==token)+obj.PriorWeight)/(numel(bad)+obj.PriorWeight*count); scores(k)=log(pg)-log(pb);
            end
            [~,index]=max(scores); value=obj.choiceAt(choices,index);
        end
        function beforeTrial(obj,study,trial) %#ok<INUSD>
        end
        function afterTrial(obj,study,trial) %#ok<INUSD>
        end
    end
    methods (Access=private)
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
        function width=bandwidth(~,values,span)
            if numel(values)<2, width=span/10; else, width=max(std(values),(max(values)-min(values))/sqrt(numel(values))); end
            width=max(width,max(span,eps)/1000);
        end
        function value=logKde(~,point,values,width)
            z=(point-values(:))/width; value=log(mean(exp(-0.5*z.^2)))-log(width)-0.5*log(2*pi);
        end
        function token=token(~,value)
            if isnumeric(value)&&isscalar(value), token="numeric:"+string(value,"%.17g"); else, token=string(jsonencode(value)); end
        end
        function value=choiceAt(~,choices,index)
            if iscell(choices), value=choices{index}; else, value=choices(index); end
        end
    end
end
