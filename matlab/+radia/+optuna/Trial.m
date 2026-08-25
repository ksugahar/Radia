classdef Trial < radia.optuna.BaseTrial
    %TRIAL One define-by-run evaluation in a table-backed Study.

    properties (SetAccess=private)
        Study
        Number (1,1) double
        State (1,1) string = "RUNNING"
        Value (1,1) double = NaN
        Values double = NaN
        Params struct = struct()
        Distributions struct = struct()
        IntermediateValues table
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        Constraints double = zeros(1,0)
        ErrorMessage (1,1) string = ""
        LastStep double = NaN
    end

    properties (Dependent, SetAccess=private)
        StartTime
        EndTime
    end

    properties (Access=private)
        StartTimeSerial (1,1) double = NaN
        EndTimeSerial (1,1) double = NaN
        RelativeParams struct = struct()
        RelativeDistributions struct = struct()
        ParameterNames struct = struct()
        ParameterDistributions struct = struct()
        FixedParameterNames string = strings(0,1)
        FixedParameterValues cell = cell(0,1)
    end

    methods (Hidden=true)
        function obj = Trial(study, number)
            obj.Study = study;
            obj.Number = number;
            obj.StartTimeSerial=now; %#ok<TNOW1> hot-path serial timestamp
            obj.IntermediateValues = ...
                radia.optuna.Trial.emptyIntermediateTable();
        end

        function markFinished(obj, state, value, endTime, message)
            obj.State = state;
            values = reshape(double(value), 1, []);
            obj.Value = values(1);
            obj.Values = values;
            obj.EndTimeSerial = double(endTime);
            obj.ErrorMessage = message;
        end

        function markRunning(obj, startTime)
            obj.State = "RUNNING";
            obj.StartTimeSerial = double(startTime);
            obj.EndTimeSerial = NaN;
            obj.ErrorMessage = "";
        end

        function restoreSnapshot(obj,snapshot)
            arguments
                obj
                snapshot (1,1) radia.optuna.FrozenTrial
            end
            obj.State=snapshot.State;
            obj.Value=snapshot.Value;
            obj.Values=snapshot.Values;
            obj.Params=snapshot.Params;
            obj.Distributions=snapshot.Distributions;
            obj.IntermediateValues=snapshot.IntermediateValues;
            obj.UserAttrs=snapshot.UserAttrs;
            obj.SystemAttrs=snapshot.SystemAttrs;
            obj.Constraints=snapshot.Constraints;
            obj.StartTimeSerial=datenum(snapshot.DatetimeStart); %#ok<DATNM>
            obj.EndTimeSerial=datenum(snapshot.DatetimeComplete); %#ok<DATNM>
            obj.ErrorMessage=snapshot.ErrorMessage;
            if isempty(snapshot.last_step())
                obj.LastStep=NaN;
            else
                obj.LastStep=snapshot.last_step();
            end
            names=string(fieldnames(snapshot.Params));
            for index=1:numel(names)
                obj.ParameterNames.(names(index))=names(index);
                if isfield(snapshot.Distributions,names(index))
                    obj.ParameterDistributions.(names(index))= ...
                        snapshot.Distributions.(names(index));
                end
            end
        end

        function setConstraints(obj, values)
            obj.Constraints = reshape(double(values), 1, []);
        end

        function value=startTimeSerial(obj)
            value=obj.StartTimeSerial;
        end

        function [names,values,distributions]=parameterRecords(obj)
            keys=string(fieldnames(obj.ParameterNames));
            keep=false(size(keys));
            for index=1:numel(keys)
                key=char(keys(index));
                keep(index)=isfield(obj.Params,key) && ...
                    isfield(obj.ParameterDistributions,key);
            end
            keys=keys(keep);
            count=numel(keys);
            names=strings(count,1);
            values=zeros(count,1);
            distributions=cell(count,1);
            for index=1:count
                key=char(keys(index));
                names(index)=string(obj.ParameterNames.(key));
                distribution=obj.ParameterDistributions.(key);
                distributions{index}=distribution;
                value=obj.Params.(key);
                if distribution.kind=="categorical"
                    token=radia.optuna.internal.DistributionCodec. ...
                        choiceToken(value);
                    tokens=radia.optuna.internal.DistributionCodec. ...
                        choiceTokens(distribution.choices);
                    choice=find(tokens==token,1);
                    if isempty(choice)
                        error("radia:optuna:CategoricalValue", ...
                            "Parameter '%s' is outside its categorical choices.", ...
                            names(index));
                    end
                    values(index)=choice;
                else
                    values(index)=double(value);
                end
            end
        end

        function setRelativeParameters(obj, searchSpace, values, source)
            if nargin < 4
                source = "tpe";
            end
            if ~iscell(values) || numel(values) ~= numel(searchSpace)
                error("radia:optuna:RelativeParameters", ...
                    "Relative parameter values must match the search space.");
            end
            names = strings(1, numel(searchSpace));
            for index = 1:numel(searchSpace)
                names(index) = searchSpace(index).name;
                key = obj.claimParameterName(searchSpace(index).name);
                obj.RelativeParams.(key) = values{index};
                obj.RelativeDistributions.(key) = searchSpace(index).distribution;
            end
            if strlength(string(source))==0
                return
            end
            attribute = matlab.lang.makeValidName( ...
                string(source) + "_relative_search_space");
            obj.setSystemAttr(attribute,names);
        end

        function setFixedParameters(obj,names,values)
            names=reshape(string(names),[],1);
            values=reshape(values,[],1);
            if numel(names)~=numel(values) || ...
                    numel(unique(names))~=numel(names)
                error("radia:optuna:FixedParameter", ...
                    "Fixed parameter names and values must align and be unique.");
            end
            obj.FixedParameterNames=names;
            obj.FixedParameterValues=values;
            obj.setSystemAttr("fixed_params",names);
        end

        function removeRelativeParameters(obj,names)
            for name=reshape(string(names),1,[])
                key=matlab.lang.makeValidName(name);
                if isfield(obj.RelativeParams,key)
                    obj.RelativeParams=rmfield(obj.RelativeParams,key);
                end
                if isfield(obj.RelativeDistributions,key)
                    obj.RelativeDistributions= ...
                        rmfield(obj.RelativeDistributions,key);
                end
            end
        end
    end

    methods (Static, Hidden=true)
        function value=emptyIntermediateTable()
            persistent initialized template
            if isempty(initialized)
                template=table('Size',[0,3], ...
                    'VariableTypes',{'double','double','datetime'}, ...
                    'VariableNames',{'Step','Value','Timestamp'});
                template.Timestamp.TimeZone="local";
                initialized=true;
            end
            value=template;
        end

    end

    methods
        function value=get.StartTime(obj)
            value=datetime(obj.StartTimeSerial,'ConvertFrom','datenum', ...
                'TimeZone','local');
        end

        function value=get.EndTime(obj)
            value=datetime(obj.EndTimeSerial,'ConvertFrom','datenum', ...
                'TimeZone','local');
        end

        function value = suggestFloat(obj, name, low, high, options)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = NaN
            end
            obj.ensureRunning();
            effectiveHigh=high;
            if isfinite(options.Step)
                count=floor((high-low)/options.Step+1e-12);
                effectiveHigh=low+count*options.Step;
                if effectiveHigh~=high
                    warning("radia:optuna:DistributionAdjusted", ...
                        "The distribution is specified by [%g, %g] and " + ...
                         "Step=%g, but the range is not divisible by Step. " + ...
                         "It will be replaced with [%g, %g].", ...
                        low,high,options.Step,low,effectiveHigh);
                end
            end
            spec = radia.optuna.internal.DistributionCodec.float( ...
                low, effectiveHigh, options.Log, options.Step);
            [exists, value, key] = obj.existingParameter(name, spec);
            if exists
                return;
            end
            [isFixed,value] = obj.fixedValue(name,spec);
            [isRelative,relativeValue] = obj.relativeValue(name, spec, key);
            if isFixed
                % Enqueued values override sampler-relative proposals.
            elseif isRelative
                value=relativeValue;
            else
                value = obj.Study.sampleFloat(obj, name, low, effectiveHigh, ...
                    struct("Log", options.Log, "Step", options.Step));
            end
            obj.Params.(key) = value;
            obj.ParameterDistributions.(key) = spec;
            obj.Distributions.(key) = struct( ...
                "name", "FloatDistribution", "low", low, "high", effectiveHigh, ...
                "log", options.Log, "step", options.Step);
            obj.Study.recordParameter(obj, name, "float", value, ...
                radia.optuna.internal.DistributionCodec.encode(spec));
        end

        function value = suggestInteger(obj, name, low, high)
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
            end
            value = obj.suggestIntegerImpl(name, low, high, 1, false);
        end

        function value = suggestCategorical(obj, name, choices)
            arguments
                obj
                name (1,1) string
                choices
            end
            obj.ensureRunning();
            spec = ...
                radia.optuna.internal.DistributionCodec.categorical(choices);
            [exists, value, key] = obj.existingParameter(name, spec);
            if exists
                return;
            end
            [isFixed,value] = obj.fixedValue(name,spec);
            [isRelative,relativeValue] = obj.relativeValue(name, spec, key);
            if isFixed
                % Enqueued values override sampler-relative proposals.
            elseif isRelative
                value=relativeValue;
            else
                value = obj.Study.sampleCategorical(obj, name, choices);
            end
            obj.Params.(key) = value;
            obj.ParameterDistributions.(key) = spec;
            obj.Distributions.(key) = struct( ...
                "name", "CategoricalDistribution", "choices", choices);
            obj.Study.recordParameter(obj, name, "categorical", value, ...
                radia.optuna.internal.DistributionCodec.encode(spec));
        end

        function values = suggestVector(obj, names, lows, highs, options)
            %SUGGESTVECTOR Register and sample a numeric joint search space.
            arguments
                obj
                names (1,:) string
                lows (1,:) double
                highs (1,:) double
                options.Log (1,:) logical = false(1, numel(names))
            end
            obj.ensureRunning();
            if isempty(names) || numel(unique(names)) ~= numel(names)
                error("radia:optuna:JointNames", ...
                    "Joint parameter names must be nonempty and unique.");
            end
            if numel(lows) ~= numel(names) || numel(highs) ~= numel(names)
                error("radia:optuna:JointShape", ...
                    "Joint names and bounds must have the same length.");
            end
            specs = repmat( ...
                radia.optuna.internal.DistributionCodec.float(0, 1, false, NaN), ...
                1, numel(names));
            keys = strings(1, numel(names));
            already = false(1, numel(keys));
            for index = 1:numel(keys)
                specs(index) = ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    lows(index), highs(index), options.Log(index), NaN);
                [already(index), ~, keys(index)] = ...
                    obj.existingParameter(names(index), specs(index));
            end
            if any(already) && ~all(already)
                error("radia:optuna:JointState", ...
                    "A joint search space cannot be partially sampled.");
            end
            if all(already)
                values = zeros(1, numel(keys));
                for index = 1:numel(keys)
                    values(index) = obj.Params.(keys(index));
                end
                return;
            end
            values = zeros(1, numel(names));
            available = false(1,numel(names));
            for index = 1:numel(names)
                [fixed,fixedValue] = obj.fixedValue(names(index),specs(index));
                [relative,relativeValue] = ...
                    obj.relativeValue(names(index), specs(index),keys(index));
                if fixed
                    values(index)=fixedValue;
                    available(index)=true;
                elseif relative
                    values(index)=relativeValue;
                    available(index)=true;
                end
            end
            if ~all(available)
                proposed = obj.Study.sampleJoint(obj, names, lows, highs, ...
                    Log=options.Log);
                values(~available)=proposed(~available);
            end
            for index = 1:numel(keys)
                obj.Params.(keys(index)) = values(index);
                obj.ParameterDistributions.(keys(index)) = specs(index);
                obj.Distributions.(keys(index)) = struct( ...
                    "name", "FloatDistribution", "low", lows(index), ...
                    "high", highs(index), "log", options.Log(index), ...
                    "step", NaN);
                obj.Study.recordParameter(obj, names(index), "float", ...
                    values(index), ...
                    radia.optuna.internal.DistributionCodec.encode(specs(index)));
            end
        end

        function values = suggest_vector(obj, names, lows, highs, options)
            arguments
                obj
                names (1,:) string
                lows (1,:) double
                highs (1,:) double
                options.Log (1,:) logical = false(1, numel(names))
            end
            values = obj.suggestVector(names, lows, highs, Log=options.Log);
        end

        function value = suggest_float(obj, name, low, high, options)
            %SUGGEST_FLOAT Optuna-compatible snake_case spelling.
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = NaN
            end
            value = obj.suggestFloat(name, low, high, ...
                Log=options.Log, Step=options.Step);
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

        function value = suggest_int(obj, name, low, high, options)
            %SUGGEST_INT Optuna-compatible integer suggestion.
            arguments
                obj
                name (1,1) string
                low (1,1) double
                high (1,1) double
                options.Step (1,1) double = 1
                options.Log (1,1) logical = false
            end
            value = obj.suggestIntegerImpl( ...
                name, low, high, options.Step, options.Log);
        end

        function value = suggest_categorical(obj, name, choices)
            %SUGGEST_CATEGORICAL Optuna-compatible snake_case spelling.
            value = obj.suggestCategorical(name, choices);
        end

        function report(obj, value, step)
            arguments
                obj
                value (1,1) double
                step (1,1) double
            end
            obj.ensureRunning();
            if numel(obj.Study.Directions) ~= 1
                error("radia:optuna:ReportMultiObjective", ...
                    "Trial.report is not supported for multi-objective studies.");
            end
            if ~isfinite(step) || step ~= floor(step) || step < 0
                error("radia:optuna:Report", ...
                    "step must be a nonnegative integer.");
            end
            row = obj.IntermediateValues.Step == step;
            if any(row)
                warning("radia:optuna:DuplicateReport", ...
                    "The reported value is ignored because step %d was already reported.", ...
                    step);
                return
            end
            obj.IntermediateValues(end+1,:) = {step, value, ...
                datetime("now", "TimeZone", "local")};
            obj.IntermediateValues = sortrows(obj.IntermediateValues, "Step");
            obj.LastStep = step;
            obj.Study.recordIntermediate(obj, value, step);
        end

        function decision = shouldPrune(obj)
            if numel(obj.Study.Directions) ~= 1
                error("radia:optuna:PrunerMultiObjective", ...
                    "Trial.shouldPrune is not supported for multi-objective studies.");
            end
            decision = obj.Study.Pruner.shouldPrune(obj.Study, obj);
        end

        function decision = should_prune(obj)
            decision = obj.shouldPrune();
        end

        function prune(obj)
            obj.ensureRunning();
            obj.Study.tell(obj, State="PRUNED");
        end

        function setUserAttr(obj, name, value)
            arguments
                obj
                name (1,1) string
                value
            end
            obj.ensureRunning();
            obj.UserAttrs.(matlab.lang.makeValidName(name)) = value;
            obj.Study.recordUserAttribute(obj, name, value);
        end

        function set_user_attr(obj, name, value)
            obj.setUserAttr(name, value);
        end

        function setSystemAttr(obj, name, value)
            arguments
                obj
                name (1,1) string
                value
            end
            obj.ensureRunning();
            obj.SystemAttrs.(matlab.lang.makeValidName(name)) = value;
            obj.Study.recordSystemAttribute(obj,name,value);
        end

        function set_system_attr(obj, name, value)
            obj.setSystemAttr(name, value);
        end

        function value = user_attrs(obj)
            value = obj.UserAttrs;
        end

        function value = system_attrs(obj)
            value = obj.SystemAttrs;
        end

        function value = params(obj)
            value = obj.Params;
        end

        function value = distributions(obj)
            value = obj.Distributions;
        end
    end

    methods (Access=private)
        function value = suggestIntegerImpl(obj, name, low, high, step, logScale)
            obj.ensureRunning();
            if ~(isfinite(low) && isfinite(high) && isfinite(step)) || ...
                    low ~= floor(low) || high ~= floor(high) || ...
                    step ~= floor(step) || low > high || step <= 0
                error("radia:optuna:IntegerDistribution", ...
                    "Integer bounds and Step must be finite integers with " + ...
                    "low <= high and Step > 0.");
            end
            if logScale && step ~= 1
                error("radia:optuna:IntegerDistribution", ...
                    "Log integer distributions require Step=1.");
            end
            effectiveHigh = low + floor((high - low) / step) * step;
            if effectiveHigh~=high
                warning("radia:optuna:DistributionAdjusted", ...
                    "The distribution is specified by [%g, %g] and " + ...
                     "Step=%g, but the range is not divisible by Step. " + ...
                     "It will be replaced with [%g, %g].", ...
                    low,high,step,low,effectiveHigh);
            end
            spec = radia.optuna.internal.DistributionCodec.integer( ...
                low, effectiveHigh, logScale, step);
            [exists, value, key] = obj.existingParameter(name, spec);
            if exists
                return
            end
            [isFixed,value] = obj.fixedValue(name,spec);
            [isRelative,relativeValue] = obj.relativeValue(name, spec, key);
            if isFixed
                % Enqueued values override sampler-relative proposals.
            elseif isRelative
                value=relativeValue;
            else
                if ismethod(obj.Study.Sampler,"sampleIntegerDetailed")
                    value=obj.Study.Sampler.sampleIntegerDetailed( ...
                        obj.Study,obj,name,low,effectiveHigh,step,logScale);
                elseif ~logScale && step == 1
                    value = obj.Study.sampleInteger( ...
                        obj, name, low, effectiveHigh);
                elseif low == effectiveHigh
                    value = low;
                else
                    value = obj.Study.sampleFloat(obj, name, low, ...
                        effectiveHigh, struct("Log", logScale, "Step", step));
                end
            end
            value = low + round((double(value) - low) / step) * step;
            value = min(max(value, low), effectiveHigh);
            obj.Params.(key) = value;
            obj.ParameterDistributions.(key) = spec;
            obj.Distributions.(key) = struct( ...
                "name", "IntDistribution", "low", low, ...
                "high", effectiveHigh, "log", logScale, "step", step);
            obj.Study.recordParameter(obj, name, "integer", value, ...
                radia.optuna.internal.DistributionCodec.encode(spec));
        end

        function [exists, value, key] = ...
                existingParameter(obj, name, distribution)
            key = obj.claimParameterName(name);
            exists = isfield(obj.Params, key);
            value = [];
            if ~exists
                return
            end
            if ~isfield(obj.ParameterDistributions, key)
                error("radia:optuna:DistributionState", ...
                    "Parameter '%s' has no registered distribution.", name);
            end
            stored = obj.ParameterDistributions.(key);
            if stored.kind~=distribution.kind
                error("radia:optuna:IncompatibleDistribution", ...
                    "Cannot set a different distribution kind for parameter '%s'.",name);
            end
            if stored.kind=="categorical"
                if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                        stored,distribution)
                    error("radia:optuna:IncompatibleDistribution", ...
                        "Categorical parameter '%s' does not support a dynamic value space.",name);
                end
            elseif stored.log~=distribution.log
                error("radia:optuna:IncompatibleDistribution", ...
                    "Cannot set a different log configuration for parameter '%s'.",name);
            elseif ~radia.optuna.internal.DistributionCodec.equivalent( ...
                    stored,distribution)
                warning("radia:optuna:InconsistentParameter", ...
                    "Inconsistent numeric distribution for parameter '%s'; " + ...
                     "the first value and distribution are used.",name);
            end
            value = obj.Params.(key);
        end

        function key = claimParameterName(obj, name)
            original = string(name);
            base=string(matlab.lang.makeValidName(original));
            if isfield(obj.ParameterNames,base)
                if string(obj.ParameterNames.(base))==original
                    key=char(base);
                    return
                end
            else
                obj.ParameterNames.(base)=original;
                key=char(base);
                return
            end
            existingKeys=string(fieldnames(obj.ParameterNames));
            for existingKey=reshape(existingKeys,1,[])
                if string(obj.ParameterNames.(existingKey))==original
                    key=char(existingKey);
                    return
                end
            end

            key=base;
            suffix=2;
            while isfield(obj.ParameterNames,key)
                key=matlab.lang.makeValidName(base+"_"+suffix);
                suffix=suffix+1;
            end
            obj.ParameterNames.(key) = original;
            key=char(key);
        end

        function [available, value] = relativeValue(obj, name, distribution,key)
            if nargin<4
                key = obj.claimParameterName(name);
            else
                key=char(key);
            end
            available = isfield(obj.RelativeParams, key) && ...
                isfield(obj.RelativeDistributions, key) && ...
                radia.optuna.internal.DistributionCodec.equivalent( ...
                obj.RelativeDistributions.(key), distribution);
            if available
                value = obj.RelativeParams.(key);
                if distribution.kind == "categorical"
                    token = ...
                        radia.optuna.internal.DistributionCodec.choiceToken(value);
                    tokens = ...
                        radia.optuna.internal.DistributionCodec.choiceTokens( ...
                        distribution.choices);
                    match = find(tokens == token, 1);
                    available = ~isempty(match);
                    if available
                        value = ...
                            radia.optuna.internal.DistributionCodec.choiceAt( ...
                            distribution.choices, match);
                    end
                end
            else
                value = NaN;
            end
        end

        function [available,value] = fixedValue(obj,name,distribution)
            index=find(obj.FixedParameterNames==string(name),1);
            if isempty(index)
                available=false;
                value=[];
                return
            end
            value=obj.FixedParameterValues{index};
            available=true;
            if ~obj.distributionContains(distribution,value)
                warning("radia:optuna:FixedParameter", ...
                    "Enqueued value for '%s' is outside the requested distribution; using the fixed value to match Optuna.", ...
                    name);
            end
        end

        function result = distributionContains(~,distribution,value)
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

        function ensureRunning(obj)
            if obj.State ~= "RUNNING"
                error("radia:optuna:TrialState", ...
                    "Trial %d is already %s.", obj.Number, obj.State);
            end
        end
    end
end
