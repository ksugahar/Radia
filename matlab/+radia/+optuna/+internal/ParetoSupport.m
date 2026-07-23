classdef ParetoSupport
    %PARETOSUPPORT Shared non-dominated ranking and observation extraction.
    methods (Static)
        function [x, values, trialNumbers] = numericObservations(study, name)
            params=study.ParamTable; rows=params.Name==string(name) & isfinite(params.ValueNumeric);
            [x,values,trialNumbers]=radia.optuna.internal.ParetoSupport.collect( ...
                study,params.TrialNumber(rows),params.ValueNumeric(rows));
        end

        function [tokens, values, trialNumbers] = categoricalObservations(study, name)
            params=study.ParamTable; rows=params.Name==string(name) & params.Kind=="categorical";
            raw=params.ValueText(rows); numeric=params.ValueNumeric(rows);
            for k=1:numel(raw)
                if isfinite(numeric(k)), raw(k)="numeric:"+string(numeric(k),"%.17g"); end
            end
            [tokens,values,trialNumbers]=radia.optuna.internal.ParetoSupport.collect( ...
                study,params.TrialNumber(rows),raw);
        end

        function [goodMask,goodWeights] = splitMOTPE(study,trialNumbers,values,nBelow)
            count=size(values,1); nBelow=max(0,min(count,nBelow));
            goodMask=false(count,1); goodWeights=zeros(count,1);
            if nBelow==0, return, end
            violations=zeros(count,1); feasible=true(count,1);
            for k=1:count
                constraints=study.constraintsForTrial(trialNumbers(k));
                if any(isnan(constraints))
                    feasible(k)=false; violations(k)=Inf;
                elseif ~isempty(constraints)
                    positive=max(constraints,0); feasible(k)=all(positive<=0);
                    violations(k)=sum(positive);
                end
            end
            feasibleIndices=find(feasible);
            takeFeasible=min(nBelow,numel(feasibleIndices));
            if takeFeasible>0
                selected=radia.optuna.internal.ParetoSupport.selectByParetoHV( ...
                    values(feasibleIndices,:),study.Directions,takeFeasible);
                goodMask(feasibleIndices(selected))=true;
            end
            remaining=nBelow-sum(goodMask);
            if remaining>0
                infeasibleIndices=find(~feasible);
                [~,order]=sortrows([violations(infeasibleIndices),trialNumbers(infeasibleIndices)],[1 2]);
                chosen=infeasibleIndices(order(1:min(remaining,numel(order))));
                goodMask(chosen)=true;
            end
            selected=find(goodMask);
            selectedFeasible=selected(feasible(selected));
            if ~isempty(selectedFeasible)
                weights=radia.optuna.internal.ParetoSupport.hypervolumeContributions( ...
                    values(selectedFeasible,:),study.Directions);
                goodWeights(selectedFeasible)=weights;
            end
            goodWeights(goodMask & goodWeights==0)=1e-12;
            goodWeights=goodWeights(goodMask);
            goodWeights=goodWeights/max(max(goodWeights),1e-12);
        end

        function [rank,crowding] = rankAndCrowding(values,directions)
            if isempty(values)
                rank=zeros(0,1); crowding=zeros(0,1); return
            end
            signs=ones(1,numel(directions)); signs(string(directions)=="maximize")=-1;
            normalized=double(values).*signs;
            count=size(normalized,1); rank=inf(count,1); remaining=true(count,1); level=1;
            while any(remaining)
                front=false(count,1); candidates=find(remaining);
                for ii=1:numel(candidates)
                    i=candidates(ii); others=candidates(candidates~=i);
                    dominated=any(all(normalized(others,:)<=normalized(i,:),2) & ...
                        any(normalized(others,:)<normalized(i,:),2));
                    front(i)=~dominated;
                end
                rank(front)=level; remaining(front)=false; level=level+1;
            end
            crowding=zeros(count,1);
            for level=unique(rank)'
                front=find(rank==level); n=numel(front);
                if n<=2, crowding(front)=Inf; continue, end
                for objective=1:size(normalized,2)
                    [ordered,order]=sort(normalized(front,objective)); indices=front(order);
                    crowding(indices([1 end]))=Inf; span=ordered(end)-ordered(1);
                    if span>0
                        for k=2:n-1
                            if isfinite(crowding(indices(k)))
                                crowding(indices(k))=crowding(indices(k)) + ...
                                    (ordered(k+1)-ordered(k-1))/span;
                            end
                        end
                    end
                end
            end
        end

        function order = preferenceOrder(values,directions)
            [rank,crowding]=radia.optuna.internal.ParetoSupport.rankAndCrowding(values,directions);
            finiteCrowding=crowding; finiteCrowding(isinf(finiteCrowding))=realmax;
            [~,order]=sortrows([rank,-finiteCrowding],[1 2]);
        end
    end

    methods (Static,Access=private)
        function [data,values,trialNumbers] = collect(study,trialNumbers,data)
            objectiveCount=numel(study.Directions); keep=false(numel(trialNumbers),1);
            values=NaN(numel(trialNumbers),objectiveCount);
            for k=1:numel(trialNumbers)
                complete=study.TrialTable.TrialNumber==trialNumbers(k) & ...
                    study.TrialTable.State=="COMPLETE";
                if ~any(complete), continue, end
                rows=study.ObjectiveTable.TrialNumber==trialNumbers(k);
                objectives=study.ObjectiveTable(rows,:);
                for j=1:height(objectives)
                    values(k,objectives.ObjectiveIndex(j))=objectives.Value(j);
                end
                keep(k)=all(isfinite(values(k,:)));
            end
            data=data(keep); values=values(keep,:); trialNumbers=trialNumbers(keep);
        end


        function selected = selectByParetoHV(values,directions,count)
            [rank,~]=radia.optuna.internal.ParetoSupport.rankAndCrowding(values,directions);
            selected=false(size(values,1),1);
            for level=reshape(unique(rank)',1,[])
                front=find(rank==level); remaining=count-sum(selected);
                if remaining<=0, break, end
                if numel(front)<=remaining
                    selected(front)=true;
                else
                    pool=front;
                    while numel(pool)>remaining
                        contributions=radia.optuna.internal.ParetoSupport.hypervolumeContributions( ...
                            values(pool,:),directions);
                        [~,remove]=min(contributions);
                        pool(remove)=[];
                    end
                    selected(pool)=true;
                    break
                end
            end
        end

        function weights = hypervolumeContributions(values,directions)
            if isempty(values), weights=zeros(0,1); return, end
            signs=ones(1,numel(directions)); signs(string(directions)=="maximize")=-1;
            points=double(values).*signs;
            worst=max(points,[],1); reference=worst+max(0.1*abs(worst),1e-12);
            total=radia.optuna.internal.ParetoSupport.hypervolume(points,reference);
            weights=zeros(size(points,1),1);
            for k=1:size(points,1)
                remaining=points; remaining(k,:)=[];
                weights(k)=max(total-radia.optuna.internal.ParetoSupport.hypervolume( ...
                    remaining,reference),1e-12);
            end
        end

        function volume = hypervolume(points,reference)
            if isempty(points), volume=0; return, end
            points=min(double(points),reference);
            points=points(all(points<reference,2),:);
            if isempty(points), volume=0; return, end
            if size(points,2)==1
                volume=max(0,reference(1)-min(points(:,1))); return
            end
            coordinates=sort(unique(points(:,1)));
            volume=0;
            for k=1:numel(coordinates)
                left=coordinates(k);
                if k<numel(coordinates), right=coordinates(k+1); else, right=reference(1); end
                if right<=left, continue, end
                active=points(points(:,1)<=left,2:end);
                volume=volume+(right-left)* ...
                    radia.optuna.internal.ParetoSupport.hypervolume(active,reference(2:end));
            end
        end
    end
end
