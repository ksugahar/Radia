classdef GroupDecomposedSearchSpace < handle
    %GROUPDECOMPOSEDSEARCHSPACE Maximal co-occurrence groups for TPE.

    properties (Access=private)
        Groups cell = cell(0,1)
        ProcessedTrials logical = false(0,1)
        SeenSignatures string = strings(0,1)
        CachedGroups cell = cell(0,1)
        CacheValid (1,1) logical = false
        CachedIncludePruned (1,1) logical = false
        CachedExcludeSingle (1,1) logical = true
        Revision (1,1) double = 0
    end

    methods
        function groups=calculate(obj,study,options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                options.IncludePruned (1,1) logical = false
                options.ExcludeSingle (1,1) logical = true
            end
            trials=study.trialData();
            parameters=study.parameterData();
            eligible=trials.State=="COMPLETE";
            if options.IncludePruned
                eligible=eligible | trials.State=="PRUNED";
            end
            trialCount=numel(trials.TrialNumber);
            if trialCount>numel(obj.ProcessedTrials)
                obj.ProcessedTrials(trialCount,1)=false;
            end
            newRows=find(eligible & ~obj.ProcessedTrials(1:trialCount));
            newNumbers=trials.TrialNumber(newRows);
            changed=false;
            for number=reshape(newNumbers,1,[])
                rows=find(parameters.TrialNumber==number)';
                names=parameters.Name(rows);
                distributions=cell(numel(rows),1);
                item=0;
                for row=rows
                    item=item+1;
                    distributions{item}= ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        parameters.Kind(row), ...
                        parameters.Distribution(row));
                end
                changed=obj.update(names,distributions) || changed;
            end
            if ~isempty(newNumbers)
                obj.ProcessedTrials(newRows)=true;
            end
            cacheMatches=obj.CacheValid && ...
                obj.CachedIncludePruned==options.IncludePruned && ...
                obj.CachedExcludeSingle==options.ExcludeSingle;
            if ~changed && cacheMatches
                groups=obj.CachedGroups;
                return
            end
            groups=obj.current(ExcludeSingle=options.ExcludeSingle);
            obj.CachedIncludePruned=options.IncludePruned;
        end

        function changed=update(obj,names,distributions)
            names=reshape(string(names),[],1);
            distributions=reshape(distributions,[],1);
            if numel(names)~=numel(distributions)
                error("radia:optuna:GroupSearchSpace", ...
                    "Names and distributions must have the same length.");
            end
            signature=strjoin(sort(names),"|");
            if any(obj.SeenSignatures==signature)
                changed=false;
                return
            end
            obj.SeenSignatures(end+1,1)=signature;
            current=radia.optuna.internal.IntersectionSearchSpace.empty();
            for index=1:numel(names)
                current(end+1,1)=struct( ...
                    "name",names(index), ...
                    "distribution",distributions{index}); %#ok<AGROW>
            end
            obj.Groups=radia.optuna.internal.GroupDecomposedSearchSpace. ...
                addDistributions(obj.Groups,current);
            obj.CacheValid=false;
            obj.Revision=obj.Revision+1;
            changed=true;
        end

        function groups=current(obj,options)
            arguments
                obj
                options.ExcludeSingle (1,1) logical = true
            end
            if obj.CacheValid && ...
                    obj.CachedExcludeSingle==options.ExcludeSingle
                groups=obj.CachedGroups;
                return
            end
            groups=obj.Groups;
            for index=1:numel(groups)
                group=groups{index};
                if options.ExcludeSingle
                    keep=false(numel(group),1);
                    for item=1:numel(group)
                        keep(item)=~radia.optuna.internal. ...
                            DistributionCodec.isSingle(group(item).distribution);
                    end
                    group=group(keep);
                end
                if ~isempty(group)
                    [~,order]=sort(string({group.name}));
                    group=group(order);
                end
                groups{index}=group;
            end
            groups=groups(~cellfun(@isempty,groups));
            obj.CachedGroups=groups;
            obj.CachedExcludeSingle=options.ExcludeSingle;
            obj.CacheValid=true;
        end

        function value=revision(obj)
            value=obj.Revision;
        end
    end

    methods (Static, Access=private)
        function next=addDistributions(groups,current)
            remaining=string({current.name});
            next=cell(0,1);
            for index=1:numel(groups)
                group=groups{index};
                names=string({group.name});
                common=ismember(names,remaining);
                if any(common)
                    next{end+1,1}=group(common); %#ok<AGROW>
                end
                if any(~common)
                    next{end+1,1}=group(~common); %#ok<AGROW>
                end
                remaining=remaining(~ismember(remaining,names));
            end
            if ~isempty(remaining)
                next{end+1,1}=current(ismember(string({current.name}), ...
                    remaining));
            end
        end
    end
end
