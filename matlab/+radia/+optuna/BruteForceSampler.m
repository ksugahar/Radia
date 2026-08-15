classdef BruteForceSampler < handle
    %BRUTEFORCESAMPLER Exhaust finite define-by-run search trees.

    properties (SetAccess=private)
        Seed (1,1) double = 0
        AvoidPrematureStop (1,1) logical = false
        Stream
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.brute-force-sampler-state.v1"
        SamplerName = "brute_force"
    end

    methods
        function obj=BruteForceSampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.AvoidPrematureStop (1,1) logical = false
            end
            obj.Seed=options.Seed;
            obj.AvoidPrematureStop=options.AvoidPrematureStop;
            obj.Stream=RandStream("mt19937ar","Seed",obj.Seed);
        end

        function beforeTrial(obj,study,~)
            obj.attach(study);
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options)
            if ~isfinite(options.Step)
                error("radia:optuna:BruteForceInfinite", ...
                    "Float distributions require a finite Step with BruteForceSampler.");
            end
            count=floor((high-low)/options.Step+1e-12);
            candidates=num2cell(low+(0:count)*options.Step);
            value=obj.selectCandidate(study,trial,name,candidates);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            candidates=num2cell(low:high);
            value=obj.selectCandidate(study,trial,name,candidates);
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            candidates=cell(1,numel(choices));
            for index=1:numel(choices)
                candidates{index}= ...
                    radia.optuna.internal.DistributionCodec.choiceAt( ...
                    choices,index);
            end
            value=obj.selectCandidate(study,trial,name,candidates);
        end

        function afterTrial(obj,study,trial)
            obj.recordState(study,trial.Number);
            paths=obj.buildPaths(study,NaN);
            excludeRunning=~obj.AvoidPrematureStop;
            [count,~]=obj.subtreeCount(paths,strings(1,0),excludeRunning);
            if count==0
                study.stopWhenOptimizing();
            end
        end
    end

    methods (Access=private)
        function value=selectCandidate(obj,study,trial,name,candidates)
            obj.attach(study);
            paths=obj.buildPaths(study,trial.Number);
            prefix=trial.Params;
            key=matlab.lang.makeValidName(name);
            counts=zeros(1,numel(candidates));
            running=false(1,numel(candidates));
            excludeRunning=~obj.AvoidPrematureStop;
            ignored=string(fieldnames(prefix));
            ignored(end+1)=key;
            for index=1:numel(candidates)
                candidatePrefix=prefix;
                candidatePrefix.(key)=candidates{index};
                selected=obj.filterPaths(paths,candidatePrefix);
                if isempty(selected)
                    counts(index)=1;
                    continue
                end
                [counts(index),running(index)]=obj.subtreeCount( ...
                    selected,ignored,excludeRunning);
            end
            if all(counts==0)
                chosen=1+floor(rand(obj.Stream)*numel(candidates));
            else
                exact=counts/sum(counts);
                flat=double(counts>0)/sum(counts>0);
                weights=0.5*exact+0.5*flat;
                if any(~running & weights>0)
                    weights(running)=0;
                end
                weights=weights/sum(weights);
                draw=rand(obj.Stream);
                chosen=find(draw<=cumsum(weights),1);
            end
            value=candidates{chosen};
            obj.recordState(study,trial.Number);
        end

        function selected=filterPaths(obj,paths,prefix)
            keep=false(1,numel(paths));
            names=string(fieldnames(prefix));
            for pathIndex=1:numel(paths)
                matches=true;
                for name=reshape(names,1,[])
                    entry=find(matlab.lang.makeValidName( ...
                        paths(pathIndex).names)==name,1);
                    if isempty(entry) || ~obj.sameValue( ...
                            paths(pathIndex).values{entry},prefix.(name))
                        matches=false;
                        break
                    end
                end
                keep(pathIndex)=matches;
            end
            selected=paths(keep);
        end

        function [count,isRunning]=subtreeCount(obj,paths,ignored,excludeRunning)
            isRunning=false;
            if isempty(paths)
                count=1;
                return
            end
            nextNames=strings(1,numel(paths));
            nextEntries=zeros(1,numel(paths));
            terminal=false(1,numel(paths));
            for index=1:numel(paths)
                keys=matlab.lang.makeValidName(paths(index).names);
                entry=find(~ismember(keys,ignored),1);
                if isempty(entry)
                    terminal(index)=true;
                else
                    nextEntries(index)=entry;
                    nextNames(index)=keys(entry);
                end
            end
            if any(terminal)
                if any(~terminal)
                    error("radia:optuna:BruteForceSearchSpace", ...
                        "A search-tree node is both terminal and expanded.");
                end
                finished=ismember([paths.state],["COMPLETE","PRUNED","FAIL"]);
                if any(finished)
                    count=0;
                    return
                end
                isRunning=any([paths.state]=="RUNNING");
                count=double(~excludeRunning || ~isRunning);
                return
            end
            if numel(unique(nextNames))~=1
                error("radia:optuna:BruteForceSearchSpace", ...
                    "Parameter order changed under the same search-tree branch.");
            end
            firstEntry=nextEntries(1);
            candidates=paths(1).candidates{firstEntry};
            for index=2:numel(paths)
                if ~obj.sameCandidates(candidates, ...
                        paths(index).candidates{nextEntries(index)})
                    error("radia:optuna:BruteForceSearchSpace", ...
                        "A finite distribution changed under the same branch.");
                end
            end
            count=0;
            childIgnored=[ignored,nextNames(1)];
            for candidate=1:numel(candidates)
                keep=false(1,numel(paths));
                for index=1:numel(paths)
                    keep(index)=obj.sameValue( ...
                        paths(index).values{nextEntries(index)}, ...
                        candidates{candidate});
                end
                if ~any(keep)
                    count=count+1;
                else
                    [childCount,~]=obj.subtreeCount( ...
                        paths(keep),childIgnored,excludeRunning);
                    count=count+childCount;
                end
            end
        end

        function paths=buildPaths(obj,study,excludedNumber)
            template=struct("number",0,"state","", ...
                "names",strings(1,0),"candidates",{{}},"values",{{}});
            states=["COMPLETE","PRUNED","RUNNING","FAIL"];
            trialRows=find(ismember(study.TrialTable.State,states))';
            paths=repmat(template,1,numel(trialRows));
            pathCount=0;
            for row=trialRows
                number=study.TrialTable.TrialNumber(row);
                if isfinite(excludedNumber) && number==excludedNumber
                    continue
                end
                paramRows=find(study.ParamTable.TrialNumber==number)';
                names=strings(1,numel(paramRows));
                candidates=cell(1,numel(paramRows));
                values=cell(1,numel(paramRows));
                for index=1:numel(paramRows)
                    paramRow=paramRows(index);
                    names(index)=study.ParamTable.Name(paramRow);
                    distribution= ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        study.ParamTable.Kind(paramRow), ...
                        study.ParamTable.Distribution(paramRow));
                    candidates{index}=obj.enumerate(distribution);
                    values{index}=obj.parameterValue(study.ParamTable(paramRow,:));
                end
                pathCount=pathCount+1;
                paths(pathCount)=struct("number",number, ...
                    "state",study.TrialTable.State(row),"names",names, ...
                    "candidates",{candidates},"values",{values});
            end
            paths=paths(1:pathCount);
        end

        function candidates=enumerate(~,distribution)
            if distribution.kind=="categorical"
                candidates=cell(1,numel(distribution.choices));
                for index=1:numel(candidates)
                    candidates{index}= ...
                        radia.optuna.internal.DistributionCodec.choiceAt( ...
                        distribution.choices,index);
                end
                return
            end
            if ~isfinite(distribution.step)
                error("radia:optuna:BruteForceInfinite", ...
                    "Float distributions require a finite Step with BruteForceSampler.");
            end
            count=floor((distribution.high-distribution.low)/ ...
                distribution.step+1e-12);
            candidates=num2cell(distribution.low+ ...
                (0:count)*distribution.step);
        end

        function value=parameterValue(~,row)
            if isfinite(row.ValueNumeric)
                value=row.ValueNumeric;
            elseif strlength(row.ValueText)>0
                value=jsondecode(row.ValueText);
            else
                value=NaN;
            end
        end

        function result=sameCandidates(obj,left,right)
            if numel(left)~=numel(right)
                result=false;
                return
            end
            result=true;
            for index=1:numel(left)
                if ~obj.sameValue(left{index},right{index})
                    result=false;
                    return
                end
            end
        end

        function result=sameValue(~,left,right)
            result=radia.optuna.internal.DistributionCodec.choiceToken(left)== ...
                radia.optuna.internal.DistributionCodec.choiceToken(right);
        end

        function attach(obj,study)
            changed=isempty(obj.AttachedStudy) || ~isequal(obj.AttachedStudy,study);
            if changed
                obj.AttachedStudy=study;
                obj.Stream=RandStream("mt19937ar","Seed",obj.Seed);
                obj.Restored=false;
            end
            if obj.Restored
                return
            end
            state=study.samplerState(obj.SamplerName,obj.StateSchema);
            if ~isempty(state)
                valid=isstruct(state) && isscalar(state) && ...
                    isfield(state,"schema") && isfield(state,"seed") && ...
                    isfield(state,"avoid_premature_stop") && ...
                    isfield(state,"random_state") && ...
                    string(state.schema)==obj.StateSchema && ...
                    double(state.seed)==obj.Seed && ...
                    logical(state.avoid_premature_stop)==obj.AvoidPrematureStop;
                if ~valid
                    error("radia:optuna:BruteForceState", ...
                        "Stored brute-force state is invalid or incompatible.");
                end
                obj.Stream.State=state.random_state;
            end
            obj.Restored=true;
        end

        function recordState(obj,study,trialNumber)
            obj.attach(study);
            state=struct("schema",obj.StateSchema,"seed",obj.Seed, ...
                "avoid_premature_stop",obj.AvoidPrematureStop, ...
                "random_state",obj.Stream.State);
            study.recordSamplerState(obj.SamplerName,obj.StateSchema, ...
                trialNumber,0,state);
        end
    end
end
