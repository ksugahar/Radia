classdef PartialFixedSampler < radia.optuna.BaseSampler
    %PARTIALFIXEDSAMPLER Fix selected parameters and delegate the remainder.

    properties (SetAccess=private)
        FixedParams (1,1) struct
        BaseSampler
    end

    methods
        function obj=PartialFixedSampler(fixedParams,baseSampler)
            arguments
                fixedParams (1,1) struct
                baseSampler
            end
            required=["sampleFloat","sampleInteger","sampleCategorical"];
            for method=required
                if ~ismethod(baseSampler,method)
                    error("radia:optuna:BaseSampler", ...
                        "BaseSampler must implement %s.",method);
                end
            end
            obj.FixedParams=fixedParams;
            obj.BaseSampler=baseSampler;
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options)
            [fixed,value]=obj.fixedValue(name);
            if fixed
                contained=isnumeric(value) && isscalar(value) && ...
                    isfinite(value) && value>=low && value<=high && ...
                    (~options.Log || value>0);
                if contained && isfinite(options.Step)
                    grid=(double(value)-low)/options.Step;
                    contained=abs(grid-round(grid))<= ...
                        1e-10*max(1,abs(grid));
                end
                obj.warnIfOutside(name,value,contained);
                return
            end
            value=obj.BaseSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            [fixed,value]=obj.fixedValue(name);
            if fixed
                contained=isnumeric(value) && isscalar(value) && ...
                    isfinite(value) && value==floor(value) && ...
                    value>=low && value<=high;
                obj.warnIfOutside(name,value,contained);
                return
            end
            value=obj.BaseSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            [fixed,value]=obj.fixedValue(name);
            if fixed
                token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                contained=ismember(token, ...
                    radia.optuna.internal.DistributionCodec.choiceTokens(choices));
                obj.warnIfOutside(name,value,contained);
                return
            end
            value=obj.BaseSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            if ismethod(obj.BaseSampler,"sampleJoint")
                values=obj.BaseSampler.sampleJoint( ...
                    study,trial,names,lows,highs,options);
            else
                values=zeros(1,numel(names));
                for index=1:numel(names)
                    values(index)=obj.BaseSampler.sampleFloat( ...
                        study,trial,names(index),lows(index),highs(index), ...
                        struct("Log",options.Log(index),"Step",NaN));
                end
            end
            for index=1:numel(names)
                [fixed,value]=obj.fixedValue(names(index));
                if fixed
                    contained=isnumeric(value) && isscalar(value) && ...
                        isfinite(value) && value>=lows(index) && ...
                        value<=highs(index) && ...
                        (~options.Log(index) || value>0);
                    obj.warnIfOutside(names(index),value,contained);
                    values(index)=value;
                end
            end
        end

        function beforeTrial(obj,study,trial)
            if ismethod(obj.BaseSampler,"beforeTrial")
                obj.BaseSampler.beforeTrial(study,trial);
            end
            trial.removeRelativeParameters(string(fieldnames(obj.FixedParams)));
        end

        function afterTrial(obj,study,trial)
            if ismethod(obj.BaseSampler,"afterTrial")
                obj.BaseSampler.afterTrial(study,trial);
            end
        end

        function searchSpace=inferRelativeSearchSpace(obj,study,trial)
            if ~ismethod(obj.BaseSampler,"inferRelativeSearchSpace")
                searchSpace=struct([]);
                return
            end
            searchSpace=obj.BaseSampler.inferRelativeSearchSpace(study,trial);
            if isempty(searchSpace)
                return
            end
            fixed=string(fieldnames(obj.FixedParams));
            searchSpace=searchSpace(~ismember([searchSpace.name],fixed));
        end
    end

    methods (Access=private)
        function [fixed,value]=fixedValue(obj,name)
            key=matlab.lang.makeValidName(name);
            fixed=isfield(obj.FixedParams,key);
            if fixed
                value=obj.FixedParams.(key);
            else
                value=[];
            end
        end

        function warnIfOutside(~,name,value,contained)
            if ~contained
                warning("radia:optuna:FixedParameter", ...
                    "Fixed parameter '%s' with value %s is outside the requested distribution.", ...
                    name,string(value));
            end
        end
    end
end
