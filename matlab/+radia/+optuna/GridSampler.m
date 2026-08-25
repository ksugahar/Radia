classdef GridSampler < radia.optuna.BaseSampler
    %GRIDSAMPLER Exhaustively evaluate a finite Cartesian product.

    properties (SetAccess=private)
        SearchSpace (1,1) struct
        Seed (1,1) double = 0
        ParameterNames string
        AllGrids cell
    end

    properties (Access=private)
        Stream
    end

    methods
        function obj=GridSampler(searchSpace,options)
            arguments
                searchSpace (1,1) struct
                options.Seed double = double.empty(1,0)
            end
            names=sort(string(fieldnames(searchSpace)));
            if isempty(names)
                error("radia:optuna:GridSearchSpace", ...
                    "SearchSpace must contain at least one parameter.");
            end
            normalized=struct();
            choices=cell(1,numel(names));
            for index=1:numel(names)
                value=searchSpace.(names(index));
                if isempty(value)
                    error("radia:optuna:GridSearchSpace", ...
                        "Grid parameter '%s' has no candidate values.",names(index));
                end
                if iscell(value)
                    choices{index}=reshape(value,1,[]);
                else
                    choices{index}=num2cell(reshape(value,1,[]));
                end
                normalized.(names(index))=choices{index};
            end
            obj.SearchSpace=normalized;
            obj.ParameterNames=names;
            obj.Seed=radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream=radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.AllGrids=obj.cartesianProduct(choices);
            order=randperm(obj.Stream,size(obj.AllGrids,1));
            obj.AllGrids=obj.AllGrids(order,:);
        end

        function beforeTrial(obj,study,trial)
            if isfield(trial.SystemAttrs,"fixed_params") || ...
                    isfield(trial.SystemAttrs,"grid_id")
                return
            end
            count=size(obj.AllGrids,1);
            if trial.Number>=0 && trial.Number<count
                gridId=trial.Number;
            else
                target=obj.unvisitedGridIds(study);
                if isempty(target)
                    warning("radia:optuna:GridExhausted", ...
                        "GridSampler is re-evaluating a configuration because the grid is exhausted.");
                    target=0:(count-1);
                end
                gridId=target(randi(obj.Stream,numel(target)));
            end
            trial.setSystemAttr("search_space",obj.SearchSpace);
            trial.setSystemAttr("grid_id",gridId);
        end

        function value=sampleFloat(obj,~,trial,name,low,high,options)
            value=obj.gridValue(trial,name);
            contained=isnumeric(value) && isscalar(value) && isfinite(value) && ...
                value>=low && value<=high && (~options.Log || value>0);
            if contained && isfinite(options.Step)
                grid=(double(value)-low)/options.Step;
                contained=abs(grid-round(grid))<=1e-10*max(1,abs(grid));
            end
            obj.warnIfOutside(name,value,contained);
        end

        function value=sampleInteger(obj,~,trial,name,low,high)
            value=obj.gridValue(trial,name);
            contained=isnumeric(value) && isscalar(value) && isfinite(value) && ...
                value==floor(value) && value>=low && value<=high;
            obj.warnIfOutside(name,value,contained);
        end

        function value=sampleCategorical(obj,~,trial,name,choices)
            value=obj.gridValue(trial,name);
            token=radia.optuna.internal.DistributionCodec.choiceToken(value);
            contained=ismember(token, ...
                radia.optuna.internal.DistributionCodec.choiceTokens(choices));
            obj.warnIfOutside(name,value,contained);
        end

        function values=sampleJoint(obj,~,trial,names,~,~,~)
            values=zeros(1,numel(names));
            for index=1:numel(names)
                value=obj.gridValue(trial,names(index));
                if ~isnumeric(value) || ~isscalar(value)
                    error("radia:optuna:GridJoint", ...
                        "Joint numeric parameter '%s' must be scalar numeric.", ...
                        names(index));
                end
                values(index)=value;
            end
        end

        function afterTrial(obj,study,~)
            if obj.isExhausted(study)
                study.stopWhenOptimizing();
            end
        end

        function result=isExhausted(obj,study)
            result=isempty(obj.unvisitedGridIds(study));
        end

        function result=is_exhausted(obj,study)
            result=obj.isExhausted(study);
        end
    end

    methods (Access=private)
        function value=gridValue(obj,trial,name)
            if ~isfield(trial.SystemAttrs,"grid_id")
                error("radia:optuna:GridEnqueue", ...
                    "All grid parameters must be supplied when enqueueTrial is used.");
            end
            index=find(obj.ParameterNames==string(name),1);
            if isempty(index)
                error("radia:optuna:GridParameter", ...
                    "Parameter '%s' is not present in SearchSpace.",name);
            end
            gridId=double(trial.SystemAttrs.grid_id);
            if gridId<0 || gridId>=size(obj.AllGrids,1)
                error("radia:optuna:GridIdentifier", ...
                    "Trial grid identifier is outside SearchSpace.");
            end
            value=obj.AllGrids{gridId+1,index};
        end

        function target=unvisitedGridIds(obj,study)
            count=size(obj.AllGrids,1);
            visited=zeros(1,0);
            running=zeros(1,0);
            rows=study.SystemAttrTable.Name=="grid_id";
            for index=find(rows)'
                number=study.SystemAttrTable.TrialNumber(index);
                trialRow=study.TrialTable.TrialNumber==number;
                if sum(trialRow)~=1
                    continue
                end
                gridId=double(jsondecode( ...
                    study.SystemAttrTable.ValueJSON(index)));
                state=study.TrialTable.State(trialRow);
                if ismember(state,["COMPLETE","PRUNED","FAIL"])
                    visited(end+1)=gridId; %#ok<AGROW>
                elseif state=="RUNNING"
                    running(end+1)=gridId; %#ok<AGROW>
                end
            end
            target=setdiff(0:(count-1),unique([visited running]),"stable");
            if isempty(target)
                target=setdiff(0:(count-1),unique(visited),"stable");
            end
        end

        function warnIfOutside(~,name,value,contained)
            if ~contained
                warning("radia:optuna:GridValue", ...
                    "Grid value %s is outside the requested distribution for '%s'.", ...
                    string(value),name);
            end
        end
    end

    methods (Static, Access=private)
        function grids=cartesianProduct(choices)
            dimensions=cellfun(@numel,choices);
            count=prod(dimensions);
            grids=cell(count,numel(choices));
            for row=0:(count-1)
                remainder=row;
                for column=numel(choices):-1:1
                    index=mod(remainder,dimensions(column))+1;
                    remainder=floor(remainder/dimensions(column));
                    grids{row+1,column}=choices{column}{index};
                end
            end
        end
    end
end
