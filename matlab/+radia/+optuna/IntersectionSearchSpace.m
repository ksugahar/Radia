classdef IntersectionSearchSpace < handle
    %INTERSECTIONSEARCHSPACE Public Optuna 4.9 intersection-space calculator.

    properties (SetAccess=private)
        IncludePruned (1,1) logical = false
    end

    methods
        function obj=IntersectionSearchSpace(options)
            arguments
                options.include_pruned (1,1) logical = false
            end
            obj.IncludePruned=options.include_pruned;
        end

        function searchSpace=calculate(obj,study,options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                options.use_cache (1,1) logical = false %#ok<INUSA>
            end
            entries=radia.optuna.internal.IntersectionSearchSpace.calculate( ...
                study,IncludePruned=obj.IncludePruned,ExcludeSingle=false);
            searchSpace=radia.optuna.IntersectionSearchSpace.toMap(entries);
        end
    end

    methods (Static)
        function searchSpace=fromTrials(trials,includePruned)
            arguments
                trials
                includePruned (1,1) logical = false
            end
            if iscell(trials)
                trials=[trials{:}];
            end
            if isempty(trials)
                searchSpace=containers.Map( ...
                    "KeyType","char","ValueType","any");
                return
            end
            if ~isa(trials,"radia.optuna.FrozenTrial")
                error("radia:optuna:SearchSpaceTrials", ...
                    "trials must contain FrozenTrial objects.");
            end
            selected=false(size(trials));
            for index=1:numel(trials)
                selected(index)=trials(index).State=="COMPLETE" || ...
                    (includePruned && trials(index).State=="PRUNED");
            end
            trials=trials(selected);
            if isempty(trials)
                searchSpace=containers.Map( ...
                    "KeyType","char","ValueType","any");
                return
            end

            names=sort(string(fieldnames(trials(1).Distributions)));
            for trialIndex=2:numel(trials)
                names=intersect(names, ...
                    string(fieldnames(trials(trialIndex).Distributions)), ...
                    "stable");
            end
            entries=radia.optuna.internal.IntersectionSearchSpace.empty();
            for name=reshape(names,1,[])
                distribution=trials(1).Distributions.(name);
                compatible=true;
                for trialIndex=2:numel(trials)
                    candidate=trials(trialIndex).Distributions.(name);
                    if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                            distribution,candidate)
                        compatible=false;
                        break
                    end
                end
                if compatible
                    entries(end+1,1)=struct( ...
                        "name",name,"distribution",distribution); %#ok<AGROW>
                end
            end
            searchSpace=radia.optuna.IntersectionSearchSpace.toMap(entries);
        end

        function searchSpace=toMap(entries)
            searchSpace=containers.Map("KeyType","char","ValueType","any");
            for index=1:numel(entries)
                searchSpace(char(entries(index).name))= ...
                    entries(index).distribution;
            end
        end
    end
end
