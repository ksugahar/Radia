classdef FrozenTrial < radia.optuna.BaseTrial
    %FROZENTRIAL Finished-trial snapshot matching the Optuna 5.0 subset.

    properties (SetAccess=private)
        Number (1,1) double = -1
        State (1,1) string = "COMPLETE"
        Value (1,1) double = NaN
        Values double = NaN
        Params struct = struct()
        Distributions struct = struct()
        IntermediateValues table = table()
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        Constraints double = zeros(1,0)
        ConstraintNames string = strings(1,0)
        ConstraintPresent (1,1) logical = false
        DatetimeStart datetime = NaT
        DatetimeComplete datetime = NaT
        Duration double = NaN
        ErrorMessage (1,1) string = ""
    end

    methods
        function obj=FrozenTrial(options)
            arguments
                options.Number (1,1) double = -1
                options.State = "COMPLETE"
                options.Value double = NaN
                options.Values double = NaN
                options.Params (1,1) struct = struct()
                options.Distributions (1,1) struct = struct()
                options.IntermediateValues table = table()
                options.UserAttrs (1,1) struct = struct()
                options.SystemAttrs (1,1) struct = struct()
                options.Constraints double = zeros(1,0)
                options.ConstraintNames string = strings(1,0)
                options.ConstraintPresent (1,1) logical = false
                options.DatetimeStart datetime = NaT
                options.DatetimeComplete datetime = NaT
                options.ErrorMessage (1,1) string = ""
            end
            state=radia.optuna.TrialState.toStorage(options.State);
            values=reshape(double(options.Values),1,[]);
            if all(isnan(values)) && ~all(isnan(options.Value))
                values=reshape(double(options.Value),1,[]);
            end
            if isempty(values), values=NaN; end
            obj.Number=options.Number;
            obj.State=state;
            obj.Values=values;
            obj.Value=values(1);
            obj.Params=options.Params;
            obj.Distributions=options.Distributions;
            obj.IntermediateValues=options.IntermediateValues;
            obj.UserAttrs=options.UserAttrs;
            obj.SystemAttrs=options.SystemAttrs;
            obj.Constraints=reshape(double(options.Constraints),1,[]);
            obj.ConstraintNames=reshape(string(options.ConstraintNames),1,[]);
            if isempty(obj.ConstraintNames) && ~isempty(obj.Constraints)
                obj.ConstraintNames=string(0:numel(obj.Constraints)-1);
            end
            if numel(obj.ConstraintNames)~=numel(obj.Constraints) || ...
                    numel(unique(obj.ConstraintNames))~=numel(obj.ConstraintNames)
                error("radia:optuna:ConstraintShape", ...
                    "Constraint names and values must align and be unique.");
            end
            obj.ConstraintPresent=options.ConstraintPresent || ...
                ~isempty(options.Constraints);
            obj.DatetimeStart=options.DatetimeStart;
            obj.DatetimeComplete=options.DatetimeComplete;
            obj.ErrorMessage=options.ErrorMessage;
            if ~isnat(obj.DatetimeStart) && ~isnat(obj.DatetimeComplete)
                obj.Duration=seconds(obj.DatetimeComplete-obj.DatetimeStart);
            end
        end

        function step=last_step(obj)
            if isempty(obj.IntermediateValues)
                step=[];
            else
                step=max(obj.IntermediateValues.Step);
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
            requested=radia.optuna.FloatDistribution( ...
                low,high,Log=options.Log,Step=options.Step);
            value=obj.suggest(name,requested);
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
            requested=radia.optuna.IntDistribution( ...
                low,high,Step=options.Step,Log=options.Log);
            value=obj.suggest(name,requested);
        end

        function value=suggestCategorical(obj,name,choices)
            requested=radia.optuna.CategoricalDistribution(choices);
            value=obj.suggest(name,requested);
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
            % Upstream FrozenTrial.report is intentionally a no-op.
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

        function setConstraint(obj,name,value)
            try
                value=double(value);
            catch
                error("radia:optuna:ConstraintType", ...
                    "Constraint '%s' must be convertible to a scalar double.",name);
            end
            if ~isscalar(value) || ~isreal(value)
                error("radia:optuna:ConstraintType", ...
                    "Constraint '%s' must be convertible to a scalar double.",name);
            end
            if isnan(value)
                error("radia:optuna:ConstraintNaN", ...
                    "Attempted to set constraint '%s', but NaN is not allowed.",name);
            end
            name=string(name);
            if any(obj.ConstraintNames==name)
                warning("radia:optuna:DuplicateConstraint", ...
                    "The constraint value is ignored because constraint '%s' is already set.",name);
                return
            end
            obj.ConstraintNames(end+1)=name;
            obj.Constraints(end+1)=value;
        end

        function set_constraint(obj,name,value)
            obj.setConstraint(name,value);
        end

        function value=constraints(obj)
            value=dictionary(reshape(obj.ConstraintNames,[],1), ...
                reshape(obj.Constraints,[],1));
        end

        function value=system_attrs(obj)
            value=obj.SystemAttrs;
        end
    end

    methods (Hidden)
        function setInternalAttribute(obj,name,value)
            % Internal storage/callback state; not an Optuna 5 public API.
            obj.SystemAttrs.(matlab.lang.makeValidName(name))=value;
        end
    end

    methods (Access=private)
        function value=suggest(obj,name,requested)
            key=matlab.lang.makeValidName(name);
            if ~isfield(obj.Params,key)
                error("radia:optuna:FrozenParameterMissing", ...
                    "The value of parameter '%s' is not present in FrozenTrial.", ...
                    name);
            end
            if ~isfield(obj.Distributions,key)
                error("radia:optuna:FrozenDistributionMissing", ...
                    "The distribution of parameter '%s' is not present in FrozenTrial.", ...
                    name);
            end
            stored=obj.Distributions.(key);
            if stored.kind~=requested.kind
                error("radia:optuna:IncompatibleDistribution", ...
                    "Cannot set a different distribution kind to the same parameter name.");
            end
            if stored.kind=="categorical"
                if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                        stored,requested)
                    error("radia:optuna:IncompatibleDistribution", ...
                        "Categorical parameters do not support a dynamic value space.");
                end
            elseif stored.log~=requested.log
                error("radia:optuna:IncompatibleDistribution", ...
                    "Cannot set a different log configuration to the same parameter name.");
            end
            value=obj.Params.(key);
            if ~obj.distributionContains(requested,value)
                warning("radia:optuna:FrozenParameter", ...
                    "The frozen value for '%s' is outside the requested distribution.", ...
                    name);
            end
        end

        function present=distributionContains(~,distribution,value)
            if distribution.kind=="categorical"
                tokens=radia.optuna.internal.DistributionCodec.choiceTokens( ...
                    distribution.choices);
                present=any(tokens== ...
                    radia.optuna.internal.DistributionCodec.choiceToken(value));
                return
            end
            present=isnumeric(value) && isscalar(value) && isfinite(value) && ...
                value>=distribution.low && value<=distribution.high;
            if present && isfinite(distribution.step)
                quotient=(double(value)-distribution.low)/distribution.step;
                present=abs(quotient-round(quotient))<=1e-10* ...
                    max(1,abs(quotient));
            end
        end
    end
end
