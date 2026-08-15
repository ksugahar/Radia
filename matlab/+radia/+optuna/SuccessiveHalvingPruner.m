classdef SuccessiveHalvingPruner < handle
    %SUCCESSIVEHALVINGPRUNER Asynchronous successive-halving (ASHA).

    properties (SetAccess=private)
        MinResource = "auto"
        ReductionFactor (1,1) double = 4
        MinEarlyStoppingRate (1,1) double = 0
        BootstrapCount (1,1) double = 0
    end

    methods
        function obj=SuccessiveHalvingPruner(options)
            arguments
                options.MinResource = "auto"
                options.ReductionFactor (1,1) double = 4
                options.MinEarlyStoppingRate (1,1) double = 0
                options.BootstrapCount (1,1) double = 0
            end
            minResource=options.MinResource;
            isAuto=isstring(minResource) && isscalar(minResource) && ...
                minResource=="auto";
            isPositiveInteger=isnumeric(minResource) && isscalar(minResource) && ...
                isfinite(minResource) && minResource>=1 && ...
                minResource==floor(minResource);
            if ~(isAuto || isPositiveInteger)
                error("radia:optuna:PrunerResource", ...
                    "MinResource must be a positive integer or 'auto'.");
            end
            if options.ReductionFactor<2 || ...
                    options.ReductionFactor~=floor(options.ReductionFactor)
                error("radia:optuna:PrunerReduction", ...
                    "ReductionFactor must be an integer of at least two.");
            end
            if options.MinEarlyStoppingRate<0 || ...
                    options.MinEarlyStoppingRate~= ...
                    floor(options.MinEarlyStoppingRate)
                error("radia:optuna:PrunerRate", ...
                    "MinEarlyStoppingRate must be a nonnegative integer.");
            end
            if options.BootstrapCount<0 || ...
                    options.BootstrapCount~=floor(options.BootstrapCount)
                error("radia:optuna:PrunerBootstrap", ...
                    "BootstrapCount must be a nonnegative integer.");
            end
            if options.BootstrapCount>0 && isAuto
                error("radia:optuna:PrunerBootstrap", ...
                    "BootstrapCount > 0 is incompatible with MinResource='auto'.");
            end
            obj.MinResource=minResource;
            obj.ReductionFactor=options.ReductionFactor;
            obj.MinEarlyStoppingRate=options.MinEarlyStoppingRate;
            obj.BootstrapCount=options.BootstrapCount;
        end

        function decision=shouldPrune(obj,study,trial)
            decision=obj.shouldPruneFiltered( ...
                study,trial,study.TrialTable.TrialNumber);
        end
    end

    methods (Hidden)
        function decision=shouldPruneFiltered(obj,study,trial,eligibleNumbers)
            decision=false;
            if isempty(trial.IntermediateValues)
                return
            end
            minResource=obj.MinResource;
            if isstring(minResource)
                minResource=obj.estimateMinResource(study,eligibleNumbers);
                if isempty(minResource)
                    return
                end
                obj.MinResource=minResource;
            end
            step=max(trial.IntermediateValues.Step);
            value=radia.optuna.internal.PrunerSupport.latestIntermediate(trial);
            rung=obj.currentRung(trial);
            while true
                promotionStep=minResource*obj.ReductionFactor^( ...
                    obj.MinEarlyStoppingRate+rung);
                if step<promotionStep
                    return
                end
                if isnan(value)
                    decision=true;
                    return
                end
                key="completed_rung_"+string(rung);
                trial.setSystemAttr(key,value);
                competing=obj.competingValues( ...
                    study,trial.Number,eligibleNumbers,key,value);
                if numel(competing)<=obj.BootstrapCount || ...
                        ~obj.isPromotable(value,competing,study.Directions(1))
                    decision=true;
                    return
                end
                rung=rung+1;
            end
        end
    end

    methods (Access=private)
        function rung=currentRung(~,trial)
            rung=0;
            while isfield(trial.SystemAttrs, ...
                    matlab.lang.makeValidName("completed_rung_"+string(rung)))
                rung=rung+1;
            end
        end

        function result=estimateMinResource(~,study,eligibleNumbers)
            complete=study.TrialTable.State=="COMPLETE" & ...
                ismember(study.TrialTable.TrialNumber,eligibleNumbers);
            rows=find(complete)';
            lastSteps=zeros(1,0);
            for row=rows
                values=study.TrialTable.IntermediateValues{row};
                if ~isempty(values)
                    lastSteps(end+1)=max(values.Step); %#ok<AGROW>
                end
            end
            if isempty(lastSteps)
                result=[];
            else
                result=max(floor(max(lastSteps)/100),1);
            end
        end

        function values=competingValues(~,study,currentNumber, ...
                eligibleNumbers,key,currentValue)
            rows=study.SystemAttrTable.Name==key & ...
                study.SystemAttrTable.TrialNumber~=currentNumber & ...
                ismember(study.SystemAttrTable.TrialNumber,eligibleNumbers);
            values=zeros(1,0);
            for index=find(rows)'
                item=jsondecode(study.SystemAttrTable.ValueJSON(index));
                values(end+1)=double(item); %#ok<AGROW>
            end
            values(end+1)=currentValue;
        end

        function result=isPromotable(obj,value,values,direction)
            values=sort(values);
            index=floor(numel(values)/obj.ReductionFactor)-1;
            if index==-1
                index=0;
            end
            if direction=="maximize"
                threshold=values(end-index);
                result=value>=threshold;
            else
                threshold=values(index+1);
                result=value<=threshold;
            end
        end
    end
end
