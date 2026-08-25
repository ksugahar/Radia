classdef FixedTrial < radia.optuna.BaseTrial
    %FIXEDTRIAL Evaluate an objective against caller-supplied parameters.

    properties (SetAccess=private)
        Number (1,1) double = 0
        Params struct = struct()
        Distributions struct = struct()
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        DatetimeStart datetime = datetime("now",TimeZone="local")
    end

    properties (Access=private)
        FixedParams struct = struct()
        ParameterNames struct = struct()
        ParameterDistributions struct = struct()
    end

    methods
        function obj=FixedTrial(params,options)
            arguments
                params (1,1) struct
                options.Number (1,1) double = 0
            end
            obj.Number=options.Number;
            obj.FixedParams=params;
            names=string(fieldnames(params));
            for index=1:numel(names)
                obj.ParameterNames.(names(index))=names(index);
            end
        end

        function value=suggestFloat(obj,name,low,high,options)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = NaN
            end
            effectiveHigh=high;
            if isfinite(options.Step)
                effectiveHigh=low+floor((high-low)/options.Step+1e-12)*options.Step;
                if effectiveHigh~=high
                    warning("radia:optuna:DistributionAdjusted", ...
                        "FloatDistribution high was adjusted from %g to %g to align with Step=%g.", ...
                        high,effectiveHigh,options.Step);
                end
            end
            spec=radia.optuna.internal.DistributionCodec.float( ...
                low,effectiveHigh,options.Log,options.Step);
            value=obj.suggest(name,spec);
            obj.Distributions.(obj.keyFor(name))=struct( ...
                "name","FloatDistribution","low",low,"high",effectiveHigh, ...
                "log",options.Log,"step",options.Step);
        end

        function value=suggestInteger(obj,name,low,high,options)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.Step (1,1) double = 1
                options.Log (1,1) logical = false
            end
            if low~=floor(low) || high~=floor(high) || ...
                    options.Step~=floor(options.Step) || options.Step<=0
                error("radia:optuna:IntegerDistribution", ...
                    "Integer bounds and Step must be integers and Step must be positive.");
            end
            if options.Log && options.Step~=1
                error("radia:optuna:IntegerDistribution", ...
                    "Log integer distributions require Step=1.");
            end
            effectiveHigh=low+floor((high-low)/options.Step)*options.Step;
            if effectiveHigh~=high
                warning("radia:optuna:DistributionAdjusted", ...
                    "IntDistribution high was adjusted from %g to %g to align with Step=%g.", ...
                    high,effectiveHigh,options.Step);
            end
            spec=radia.optuna.internal.DistributionCodec.integer( ...
                low,effectiveHigh,options.Log,options.Step);
            value=obj.suggest(name,spec);
            obj.Distributions.(obj.keyFor(name))=struct( ...
                "name","IntDistribution","low",low,"high",effectiveHigh, ...
                "log",options.Log,"step",options.Step);
        end

        function value=suggestCategorical(obj,name,choices)
            spec=radia.optuna.internal.DistributionCodec.categorical(choices);
            value=obj.suggest(name,spec);
            tokens=radia.optuna.internal.DistributionCodec.choiceTokens(choices);
            token=radia.optuna.internal.DistributionCodec.choiceToken(value);
            index=find(tokens==token,1);
            if isempty(index)
                error("radia:optuna:FixedParameter", ...
                    "The fixed value for '%s' is not in the categorical choices.",name);
            end
            value=radia.optuna.internal.DistributionCodec.choiceAt(choices,index);
            key=obj.keyFor(name);
            obj.Params.(key)=value;
            obj.Distributions.(key)=struct( ...
                "name","CategoricalDistribution","choices",{choices});
        end

        function value=suggest_float(obj,name,low,high,options)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.step (1,1) double = NaN
                options.log (1,1) logical = false
            end
            value=obj.suggestFloat(name,low,high, ...
                Step=options.step,Log=options.log);
        end

        function value=suggest_int(obj,name,low,high,options)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.step (1,1) double = 1
                options.log (1,1) logical = false
            end
            value=obj.suggestInteger(name,low,high, ...
                Step=options.step,Log=options.log);
        end

        function value=suggest_categorical(obj,name,choices)
            value=obj.suggestCategorical(name,choices);
        end

        function value=suggest_uniform(obj,name,low,high)
            value=obj.suggestFloat(name,low,high);
        end

        function value=suggest_loguniform(obj,name,low,high)
            value=obj.suggestFloat(name,low,high,Log=true);
        end

        function value=suggest_discrete_uniform(obj,name,low,high,q)
            value=obj.suggestFloat(name,low,high,Step=q);
        end

        function report(~,~,~)
            % FixedTrial deliberately ignores intermediate reports.
        end

        function decision=shouldPrune(~)
            decision=false;
        end

        function decision=should_prune(obj)
            decision=obj.shouldPrune();
        end

        function setUserAttr(obj,name,value)
            obj.UserAttrs.(matlab.lang.makeValidName(name))=value;
        end

        function set_user_attr(obj,name,value)
            obj.setUserAttr(name,value);
        end

        function setSystemAttr(obj,name,value)
            obj.SystemAttrs.(matlab.lang.makeValidName(name))=value;
        end

        function set_system_attr(obj,name,value)
            obj.setSystemAttr(name,value);
        end

        function value=params(obj)
            value=obj.Params;
        end

        function value=distributions(obj)
            value=obj.Distributions;
        end

        function value=user_attrs(obj)
            value=obj.UserAttrs;
        end

        function value=system_attrs(obj)
            value=obj.SystemAttrs;
        end
    end

    methods (Access=private)
        function value=suggest(obj,name,distribution)
            key=obj.keyFor(name);
            if ~isfield(obj.FixedParams,key)
                error("radia:optuna:FixedParameterMissing", ...
                    "The value of parameter '%s' was not supplied to FixedTrial.",name);
            end
            value=obj.FixedParams.(key);
            if ~obj.distributionContains(distribution,value)
                warning("radia:optuna:FixedParameter", ...
                    "The fixed value for '%s' is outside the requested distribution.",name);
            end
            if isfield(obj.ParameterDistributions,key)
                stored=obj.ParameterDistributions.(key);
                if stored.kind~=distribution.kind
                    error("radia:optuna:IncompatibleDistribution", ...
                        "Cannot set a different distribution kind to the same parameter name.");
                end
                if stored.kind=="categorical" && ...
                        ~radia.optuna.internal.DistributionCodec.equivalent( ...
                        stored,distribution)
                    error("radia:optuna:IncompatibleDistribution", ...
                        "Categorical parameters do not support a dynamic value space.");
                end
                if stored.kind~="categorical" && stored.log~=distribution.log
                    error("radia:optuna:IncompatibleDistribution", ...
                        "Cannot set a different log configuration to the same parameter name.");
                end
            end
            obj.ParameterDistributions.(key)=distribution;
            obj.Params.(key)=value;
        end

        function key=keyFor(obj,name)
            original=string(name);
            keys=string(fieldnames(obj.ParameterNames));
            for candidate=reshape(keys,1,[])
                if string(obj.ParameterNames.(candidate))==original
                    key=char(candidate);
                    return
                end
            end
            key=char(matlab.lang.makeValidName(original));
        end

        function result=distributionContains(~,distribution,value)
            if distribution.kind=="categorical"
                token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                result=ismember(token, ...
                    radia.optuna.internal.DistributionCodec.choiceTokens( ...
                    distribution.choices));
                return
            end
            result=isnumeric(value) && isscalar(value) && isfinite(value) && ...
                value>=distribution.low && value<=distribution.high && ...
                (~distribution.log || value>0);
            if result && isfinite(distribution.step)
                grid=(double(value)-distribution.low)/distribution.step;
                result=abs(grid-round(grid))<=1e-10*max(1,abs(grid));
            end
            if result && distribution.kind=="integer"
                result=double(value)==round(double(value));
            end
        end
    end
end
