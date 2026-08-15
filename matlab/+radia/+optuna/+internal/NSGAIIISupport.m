classdef NSGAIIISupport
    %NSGAIIISUPPORT Deterministic Optuna 4.9 NSGA-III transforms.

    methods (Static)
        function points=defaultReferencePoints(nObjectives,dividingParameter)
            arguments
                nObjectives (1,1) double {mustBeInteger,mustBePositive}
                dividingParameter (1,1) double ...
                    {mustBeInteger,mustBePositive} = 3
            end
            points=radia.optuna.internal.NSGAIIISupport. ...
                weakCompositions(dividingParameter,nObjectives);
        end

        function normalized=normalizeObjectives(values)
            % Match _filter_inf + _normalize_objective_values in Optuna.
            working=double(values);
            if ~ismatrix(working) || isempty(working) || any(isnan(working),"all")
                error("radia:optuna:NSGAIIIObservations", ...
                    "NSGA-III objective values must form a nonempty matrix without NaN.");
            end
            for objective=1:size(working,2)
                finiteValues=working(isfinite(working(:,objective)),objective);
                if isempty(finiteValues)
                    error("radia:optuna:NSGAIIIObservations", ...
                        "Each objective must contain at least one finite value.");
                end
                minimum=min(finiteValues);
                maximum=max(finiteValues);
                margin=3*(maximum-minimum);
                working(working(:,objective)==-Inf,objective)=minimum-margin;
                working(working(:,objective)==Inf,objective)=maximum+margin;
            end
            translated=working-min(working,[],1);
            dimensions=size(translated,2);
            extreme=zeros(dimensions,dimensions);
            for axis=1:dimensions
                weights=repmat(1e6,1,dimensions);
                weights(axis)=1;
                [~,row]=min(max(translated.*weights,[],2));
                extreme(axis,:)=translated(row,:);
            end
            usePlane=all(isfinite(extreme),"all") && ...
                rank(extreme)==dimensions;
            if usePlane
                inverseIntercept=extreme\ones(dimensions,1);
                inverseIntercept(~isfinite(inverseIntercept))=1;
            else
                intercept=max(translated,[],1);
                intercept(intercept==0)=1;
                inverseIntercept=reshape(1./intercept,[],1);
                inverseIntercept(~isfinite(inverseIntercept))=1;
            end
            normalized=translated.*reshape(inverseIntercept,1,[]);
        end

        function [associations,distances]=associate(points,references)
            points=double(points);
            references=double(references);
            if size(points,2)~=size(references,2) || ...
                    any(~isfinite(points),"all") || ...
                    any(~isfinite(references),"all") || ...
                    any(vecnorm(references,2,2)==0)
                error("radia:optuna:NSGAIIIReferencePoints", ...
                    "Finite nonzero reference rows must match the objective dimension.");
            end
            references=references./vecnorm(references,2,2);
            allDistances=zeros(size(points,1),size(references,1));
            for reference=1:size(references,1)
                direction=references(reference,:);
                projection=(points*direction')*direction;
                allDistances(:,reference)=vecnorm(points-projection,2,2);
            end
            [distances,associations]=min(allDistances,[],2);
        end
    end

    methods (Static,Access=private)
        function points=weakCompositions(total,width)
            if width==1
                points=total;
                return
            end
            rows=nchoosek(total+width-1,width-1);
            points=zeros(rows,width);
            cursor=1;
            for first=total:-1:0
                tail=radia.optuna.internal.NSGAIIISupport. ...
                    weakCompositions(total-first,width-1);
                count=size(tail,1);
                points(cursor:cursor+count-1,:)=[ ...
                    repmat(first,count,1),tail];
                cursor=cursor+count;
            end
        end
    end
end
