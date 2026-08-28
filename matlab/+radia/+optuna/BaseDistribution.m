classdef (Abstract) BaseDistribution
    %BASEDISTRIBUTION Public base class for Optuna 4.9 distributions.

    properties
        name (1,1) string = "BaseDistribution"
        kind (1,1) string = ""
        low (1,1) double = NaN
        high (1,1) double = NaN
        log (1,1) logical = false
        step (1,1) double = NaN
        choices = cell(0,1)
    end

    methods
        function result=single(obj)
            if obj.kind=="categorical"
                result=isscalar(obj.choices);
            else
                result=obj.low==obj.high;
            end
        end

        function value=to_external_repr(obj,value)
            if obj.kind=="categorical"
                index=double(value)+1;
                if ~isscalar(index) || ~isfinite(index) || ...
                        index~=fix(index) || index<1 || index>numel(obj.choices)
                    error("radia:optuna:DistributionValue", ...
                        "Categorical internal representation is out of range.");
                end
                value=radia.optuna.internal.DistributionCodec.choiceAt( ...
                    obj.choices,index);
            elseif obj.kind=="integer"
                value=fix(double(value));
            end
        end

        function value=to_internal_repr(obj,value)
            if obj.kind=="categorical"
                token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                tokens=radia.optuna.internal.DistributionCodec.choiceTokens( ...
                    obj.choices);
                index=find(tokens==token,1);
                if isempty(index)
                    error("radia:optuna:DistributionValue", ...
                        "Value is not contained in the categorical choices.");
                end
                value=index-1;
            elseif ischar(value) || (isstring(value) && isscalar(value))
                value=str2double(string(value));
            else
                value=double(value);
            end
        end

        function spec=toStruct(obj)
            spec=struct( ...
                "name",obj.name, ...
                "kind",obj.kind, ...
                "low",obj.low, ...
                "high",obj.high, ...
                "log",obj.log, ...
                "step",obj.step, ...
                "choices",{obj.choices});
        end
    end
end
