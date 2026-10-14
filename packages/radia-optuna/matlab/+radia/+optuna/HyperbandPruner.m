classdef HyperbandPruner < radia.optuna.BasePruner
    %HYPERBANDPRUNER Optuna-compatible asynchronous Hyperband brackets.

    properties (SetAccess=private)
        MinResource (1,1) double = 1
        MaxResource = "auto"
        ReductionFactor (1,1) double = 3
        BootstrapCount (1,1) double = 0
    end

    properties (Access=private)
        Pruners cell = cell(1,0)
        AllocationBudgets double = zeros(1,0)
        TotalAllocationBudget (1,1) double = 0
        NBrackets (1,1) double = 0
    end

    methods
        function obj=HyperbandPruner(options)
            arguments
                options.MinResource (1,1) double = 1
                options.MaxResource = "auto"
                options.ReductionFactor (1,1) double = 3
                options.BootstrapCount (1,1) double = 0
            end
            if options.MinResource<1 || ...
                    options.MinResource~=floor(options.MinResource)
                error("radia:optuna:PrunerResource", ...
                    "MinResource must be a positive integer.");
            end
            maxResource=options.MaxResource;
            isAuto=isstring(maxResource) && isscalar(maxResource) && ...
                maxResource=="auto";
            validMax=isnumeric(maxResource) && isscalar(maxResource) && ...
                isfinite(maxResource) && maxResource>=options.MinResource && ...
                maxResource==floor(maxResource);
            if ~(isAuto || validMax)
                error("radia:optuna:PrunerResource", ...
                    "MaxResource must be an integer >= MinResource or 'auto'.");
            end
            if options.ReductionFactor<2 || ...
                    options.ReductionFactor~=floor(options.ReductionFactor)
                error("radia:optuna:PrunerReduction", ...
                    "ReductionFactor must be an integer of at least two.");
            end
            if options.BootstrapCount<0 || ...
                    options.BootstrapCount~=floor(options.BootstrapCount)
                error("radia:optuna:PrunerBootstrap", ...
                    "BootstrapCount must be a nonnegative integer.");
            end
            if options.BootstrapCount>0 && isAuto
                error("radia:optuna:PrunerBootstrap", ...
                    "BootstrapCount > 0 is incompatible with MaxResource='auto'.");
            end
            obj.MinResource=options.MinResource;
            obj.MaxResource=maxResource;
            obj.ReductionFactor=options.ReductionFactor;
            obj.BootstrapCount=options.BootstrapCount;
        end

        function decision=shouldPrune(obj,study,trial)
            decision=false;
            obj.tryInitialize(study);
            if isempty(obj.Pruners)
                return
            end
            bracket=obj.bracketId(study,trial.Number);
            numbers=study.TrialTable.TrialNumber;
            keep=false(size(numbers));
            for index=1:numel(numbers)
                keep(index)=obj.bracketId(study,numbers(index))==bracket;
            end
            decision=obj.Pruners{bracket+1}.shouldPruneFiltered( ...
                study,trial,numbers(keep));
        end
    end

    methods (Hidden)
        function bracket=bracketId(obj,study,trialNumber)
            if isempty(obj.Pruners)
                bracket=0;
                return
            end
            checksum=obj.crc32(study.Name+"_"+string(trialNumber));
            allocation=mod(double(checksum),obj.TotalAllocationBudget);
            for index=1:obj.NBrackets
                allocation=allocation-obj.AllocationBudgets(index);
                if allocation<0
                    bracket=index-1;
                    return
                end
            end
            error("radia:optuna:HyperbandBracket", ...
                "Could not assign a Hyperband bracket.");
        end
    end

    methods (Access=private)
        function tryInitialize(obj,study)
            if ~isempty(obj.Pruners)
                return
            end
            maxResource=obj.MaxResource;
            if isstring(maxResource)
                complete=study.TrialTable.State=="COMPLETE";
                lastSteps=zeros(1,0);
                for row=find(complete)'
                    values=study.TrialTable.IntermediateValues{row};
                    if ~isempty(values)
                        lastSteps(end+1)=max(values.Step); %#ok<AGROW>
                    end
                end
                if isempty(lastSteps)
                    return
                end
                maxResource=max(lastSteps)+1;
                obj.MaxResource=maxResource;
            end
            obj.NBrackets=floor(log(maxResource/obj.MinResource)/ ...
                log(obj.ReductionFactor))+1;
            for bracket=0:(obj.NBrackets-1)
                s=obj.NBrackets-1-bracket;
                budget=ceil(obj.NBrackets*obj.ReductionFactor^s/(s+1));
                obj.AllocationBudgets(end+1)=budget;
                obj.Pruners{end+1}=radia.optuna.SuccessiveHalvingPruner( ...
                    MinResource=obj.MinResource, ...
                    ReductionFactor=obj.ReductionFactor, ...
                    MinEarlyStoppingRate=bracket, ...
                    BootstrapCount=obj.BootstrapCount);
            end
            obj.TotalAllocationBudget=sum(obj.AllocationBudgets);
        end
    end

    methods (Static, Access=private)
        function checksum=crc32(value)
            bytes=unicode2native(char(value),"UTF-8");
            checksum=uint32(hex2dec("FFFFFFFF"));
            polynomial=uint32(hex2dec("EDB88320"));
            for byte=reshape(bytes,1,[])
                checksum=bitxor(checksum,uint32(byte));
                for bit=1:8
                    if bitand(checksum,uint32(1))~=0
                        checksum=bitxor(bitshift(checksum,-1),polynomial);
                    else
                        checksum=bitshift(checksum,-1);
                    end
                end
            end
            checksum=bitcmp(checksum);
        end
    end
end
