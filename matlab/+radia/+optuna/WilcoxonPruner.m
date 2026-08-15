classdef WilcoxonPruner < handle
    %WILCOXONPRUNER Paired signed-rank pruning across problem instances.

    properties (SetAccess=private)
        PThreshold (1,1) double = 0.1
        NStartupSteps (1,1) double = 2
    end

    methods
        function obj=WilcoxonPruner(options)
            arguments
                options.PThreshold (1,1) double = 0.1
                options.NStartupSteps (1,1) double = 2
            end
            if options.PThreshold<0 || options.PThreshold>1 || ...
                    ~isfinite(options.PThreshold)
                error("radia:optuna:WilcoxonThreshold", ...
                    "PThreshold must be between zero and one inclusive.");
            end
            if options.NStartupSteps<0 || ...
                    options.NStartupSteps~=floor(options.NStartupSteps)
                error("radia:optuna:WilcoxonStartup", ...
                    "NStartupSteps must be a nonnegative integer.");
            end
            obj.PThreshold=options.PThreshold;
            obj.NStartupSteps=options.NStartupSteps;
        end

        function decision=shouldPrune(obj,study,trial)
            decision=false;
            if isempty(trial.IntermediateValues)
                return
            end
            current=sortrows(trial.IntermediateValues,"Step");
            if any(~isfinite(current.Value))
                warning("radia:optuna:WilcoxonNonfinite", ...
                    "Current trial has nonfinite intermediate values and will not be pruned.");
                return
            end
            try
                bestRow=study.bestTrial();
            catch exception
                if exception.identifier=="radia:optuna:NoCompletedTrials" || ...
                        exception.identifier=="radia:optuna:NoFeasibleTrial"
                    return
                end
                rethrow(exception)
            end
            best=study.freezeTrial(bestRow.TrialNumber);
            if isempty(best.IntermediateValues)
                warning("radia:optuna:WilcoxonNoIntermediate", ...
                    "The best trial has no intermediate values.");
                return
            end
            reference=sortrows(best.IntermediateValues,"Step");
            if any(~isfinite(reference.Value))
                warning("radia:optuna:WilcoxonNonfinite", ...
                    "The best trial has nonfinite intermediate values.");
                return
            end
            [~,currentIndex,referenceIndex]=intersect( ...
                current.Step,reference.Step,"stable");
            differences=current.Value(currentIndex)- ...
                reference.Value(referenceIndex);
            if numel(differences)<max(2,obj.NStartupSteps)
                return
            end
            if study.Directions(1)=="maximize"
                alternative="less";
                safety=mean(reference.Value)<=mean(current.Value);
            else
                alternative="greater";
                safety=mean(reference.Value)>=mean(current.Value);
            end
            probability=obj.signedRankPValue(differences,alternative);
            if probability<obj.PThreshold && safety
                return
            end
            decision=probability<obj.PThreshold;
        end
    end

    methods (Access=private)
        function probability=signedRankPValue(obj,differences,alternative)
            ranks=obj.averageRanks(abs(differences));
            nonzero=differences~=0;
            weights=round(4*ranks(nonzero));
            observed=sum(weights(differences(nonzero)>0));
            count=numel(weights);
            if count<=50
                distribution=1;
                for weight=reshape(weights,1,[])
                    shifted=zeros(1,numel(distribution)+weight);
                    shifted(1:numel(distribution))=distribution;
                    shifted((weight+1):(weight+numel(distribution)))= ...
                        shifted((weight+1):(weight+numel(distribution)))+ ...
                        distribution;
                    distribution=shifted;
                end
                support=0:(numel(distribution)-1);
                if alternative=="greater"
                    probability=sum(distribution(support>=observed))/2^count;
                else
                    probability=sum(distribution(support<=observed))/2^count;
                end
                return
            end
            positive=sum(ranks(differences>0))+0.5*sum(ranks(~nonzero));
            center=0.5*sum(ranks);
            deviation=sqrt(0.25*sum(ranks(nonzero).^2));
            if deviation==0
                probability=1;
                return
            end
            z=(positive-center)/deviation;
            if alternative=="greater"
                probability=0.5*erfc(z/sqrt(2));
            else
                probability=0.5*erfc(-z/sqrt(2));
            end
        end

        function ranks=averageRanks(~,values)
            [ordered,order]=sort(values);
            ranked=zeros(size(values));
            first=1;
            while first<=numel(ordered)
                last=first;
                while last<numel(ordered) && ordered(last+1)==ordered(first)
                    last=last+1;
                end
                ranked(order(first:last))=(first+last)/2;
                first=last+1;
            end
            ranks=ranked;
        end
    end
end
