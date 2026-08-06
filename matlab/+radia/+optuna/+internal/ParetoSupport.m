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
            if ~isempty(study.ConstraintTable)
                for k=1:count
                    constraints=study.constraintsForTrial(trialNumbers(k));
                    if any(isnan(constraints))
                        feasible(k)=false; violations(k)=Inf;
                    elseif ~isempty(constraints)
                        positive=max(constraints,0); feasible(k)=all(positive<=0);
                        violations(k)=sum(positive);
                    end
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
            values=double(values); directions=string(directions);
            if ~ismatrix(values) || size(values,2)~=numel(directions) || ...
                    any(~isfinite(values),"all")
                error("radia:optuna:InvalidParetoValues", ...
                    "Values must be a finite matrix with one column per direction.");
            end
            if any(directions~="minimize" & directions~="maximize")
                error("radia:optuna:InvalidDirection", ...
                    "Directions must be minimize or maximize.");
            end
            signs=ones(1,numel(directions)); signs(directions=="maximize")=-1;
            if radia.optuna.internal.NativeKernels.has( ...
                    "optuna.pareto.rank_crowding")
                [rank,crowding]=radia.internal.callMex( ...
                    "optuna.pareto.rank_crowding",values,signs);
                return
            end
            normalized=values.*signs;
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
            objectiveCount = numel(study.Directions);
            count = numel(trialNumbers);
            values = NaN(count, objectiveCount);
            [knownTrials, trialRows] = ismember( ...
                trialNumbers, study.TrialTable.TrialNumber);
            complete = false(count, 1);
            complete(knownTrials) = ...
                study.TrialTable.State(trialRows(knownTrials)) == "COMPLETE";
            objectives = study.ObjectiveTable;
            for objective = 1:objectiveCount
                rows = objectives.ObjectiveIndex == objective;
                numbers = objectives.TrialNumber(rows);
                objectiveValues = objectives.Value(rows);
                [present, locations] = ismember(trialNumbers, numbers);
                values(present, objective) = objectiveValues(locations(present));
            end
            keep = complete & all(isfinite(values), 2);
            data = data(keep);
            values = values(keep,:);
            trialNumbers = trialNumbers(keep);
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
                    chosen = ...
                        radia.optuna.internal.ParetoSupport.greedyHVSubset( ...
                        values(front,:), directions, remaining);
                    selected(front(chosen)) = true;
                    break
                end
            end
        end

        function selected = greedyHVSubset(values, directions, count)
            % Match Optuna HSSP's forward greedy hypervolume selection.
            signs = ones(1, numel(directions));
            signs(string(directions) == "maximize") = -1;
            points = double(values) .* signs;
            [points, order] = sortrows(points);
            worst = max(points, [], 1);
            reference = worst + max(0.1 * abs(worst), 1e-12);
            if size(points, 2) == 2
                selected = ...
                    radia.optuna.internal.ParetoSupport.greedyHVSubset2D( ...
                    points, order, count, reference);
                return
            end
            chosen = false(size(points, 1), 1);
            currentVolume = 0;
            for slot = 1:min(count, size(points, 1))
                contributions = -Inf(size(points, 1), 1);
                for candidate = reshape(find(~chosen), 1, [])
                    included = chosen;
                    included(candidate) = true;
                    contributions(candidate) = ...
                        radia.optuna.internal.ParetoSupport.hypervolume( ...
                        points(included,:), reference) - currentVolume;
                end
                [~, best] = max(contributions);
                chosen(best) = true;
                currentVolume = ...
                    radia.optuna.internal.ParetoSupport.hypervolume( ...
                    points(chosen,:), reference);
            end
            selected = sort(order(chosen));
        end

        function selected = greedyHVSubset2D( ...
                points, originalOrder, count, reference)
            % O(subset size * population) rectangle update for two objectives.
            remainingIndices = (1:size(points, 1)).';
            rectangleDiagonals = repmat(reference, size(points, 1), 1);
            selectedSorted = zeros(min(count, size(points, 1)), 1);
            for slot = 1:numel(selectedSorted)
                contributions = prod(rectangleDiagonals - points, 2);
                [~, best] = max(contributions);
                selectedSorted(slot) = remainingIndices(best);
                selectedPoint = points(best, :);

                keep = true(size(points, 1), 1);
                keep(best) = false;
                remainingIndices = remainingIndices(keep);
                rectangleDiagonals = rectangleDiagonals(keep, :);
                points = points(keep, :);
                if best > 1
                    rectangleDiagonals(1:best-1, 1) = min( ...
                        selectedPoint(1), rectangleDiagonals(1:best-1, 1));
                end
                if best <= size(rectangleDiagonals, 1)
                    rectangleDiagonals(best:end, 2) = min( ...
                        selectedPoint(2), rectangleDiagonals(best:end, 2));
                end
            end
            selected = sort(originalOrder(selectedSorted));
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
