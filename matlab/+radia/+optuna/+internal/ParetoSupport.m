classdef ParetoSupport
    %PARETOSUPPORT Shared non-dominated ranking and observation extraction.
    methods (Static)
        function [x, values] = numericObservations(study, name)
            params=study.ParamTable; rows=params.Name==string(name) & isfinite(params.ValueNumeric);
            [x,values]=radia.optuna.internal.ParetoSupport.collect( ...
                study,params.TrialNumber(rows),params.ValueNumeric(rows));
        end

        function [tokens, values] = categoricalObservations(study, name)
            params=study.ParamTable; rows=params.Name==string(name) & params.Kind=="categorical";
            raw=params.ValueText(rows); numeric=params.ValueNumeric(rows);
            for k=1:numel(raw)
                if isfinite(numeric(k)), raw(k)="numeric:"+string(numeric(k),"%.17g"); end
            end
            [tokens,values]=radia.optuna.internal.ParetoSupport.collect( ...
                study,params.TrialNumber(rows),raw);
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
        function [data,values] = collect(study,trialNumbers,data)
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
            data=data(keep); values=values(keep,:);
        end
    end
end
