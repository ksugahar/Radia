classdef Trial < handle
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
        StartTime datetime
        EndTime datetime = NaT
        ErrorMessage (1,1) string = ""
        LastStep double = NaN
    end

    properties (Access=private)
        RelativeParams struct = struct()
        RelativeDistributions struct = struct()
    end

    methods (Hidden=true)
        function obj = Trial(study, number)
            obj.Study = study;
            obj.Number = number;
            obj.StartTime = datetime("now", "TimeZone", "local");
            obj.IntermediateValues = table('Size', [0, 3], ...
                'VariableTypes', {'double', 'double', 'datetime'}, ...
                'VariableNames', {'Step', 'Value', 'Timestamp'});
            obj.IntermediateValues.Timestamp.TimeZone = "local";
        end

        function markFinished(obj, state, value, endTime, message)
            obj.State = state;
            values = reshape(double(value), 1, []);
            obj.Value = values(1);
            obj.Values = values;
            obj.EndTime = endTime;
            obj.ErrorMessage = message;
        end

        function setConstraints(obj, values)
            obj.Constraints = reshape(double(values), 1, []);
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
                key = matlab.lang.makeValidName(searchSpace(index).name);
                obj.RelativeParams.(key) = values{index};
                obj.RelativeDistributions.(key) = searchSpace(index).distribution;
            end
            attribute = matlab.lang.makeValidName( ...
                string(source) + "_relative_search_space");
            obj.SystemAttrs.(attribute) = names;
        end
    end

    methods
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
            key = matlab.lang.makeValidName(name);
            if isfield(obj.Params, key)
                value = obj.Params.(key);
                return;
            end
            spec = radia.optuna.internal.DistributionCodec.float( ...
                low, high, options.Log, options.Step);
            [isRelative, value] = obj.relativeValue(name, spec);
            if ~isRelative
                value = obj.Study.sampleFloat(obj, name, low, high, ...
                    struct("Log", options.Log, "Step", options.Step));
            end
            obj.Params.(key) = value;
            obj.Distributions.(key) = struct( ...
                "name", "FloatDistribution", "low", low, "high", high, ...
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
            obj.ensureRunning();
            key = matlab.lang.makeValidName(name);
            if isfield(obj.Params, key)
                value = obj.Params.(key);
                return;
            end
            spec = radia.optuna.internal.DistributionCodec.integer( ...
                low, high, false, 1);
            [isRelative, value] = obj.relativeValue(name, spec);
            if ~isRelative
                value = obj.Study.sampleInteger(obj, name, low, high);
            end
            obj.Params.(key) = value;
            obj.Distributions.(key) = struct( ...
                "name", "IntDistribution", "low", low, "high", high, ...
                "log", false, "step", 1);
            obj.Study.recordParameter(obj, name, "integer", value, ...
                radia.optuna.internal.DistributionCodec.encode(spec));
        end

        function value = suggestCategorical(obj, name, choices)
            arguments
                obj
                name (1,1) string
                choices
            end
            obj.ensureRunning();
            key = matlab.lang.makeValidName(name);
            if isfield(obj.Params, key)
                value = obj.Params.(key);
                return;
            end
            spec = ...
                radia.optuna.internal.DistributionCodec.categorical(choices);
            [isRelative, value] = obj.relativeValue(name, spec);
            if ~isRelative
                value = obj.Study.sampleCategorical(obj, name, choices);
            end
            obj.Params.(key) = value;
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
            keys = matlab.lang.makeValidName(names);
            already = false(1, numel(keys));
            for index = 1:numel(keys)
                already(index) = isfield(obj.Params, keys(index));
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
            specs = repmat( ...
                radia.optuna.internal.DistributionCodec.float(0, 1, false, NaN), ...
                1, numel(names));
            values = zeros(1, numel(names));
            allRelative = true;
            for index = 1:numel(names)
                specs(index) = ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    lows(index), highs(index), options.Log(index), NaN);
                [available, values(index)] = ...
                    obj.relativeValue(names(index), specs(index));
                allRelative = allRelative && available;
            end
            if ~allRelative
                values = obj.Study.sampleJoint(obj, names, lows, highs, ...
                    Log=options.Log);
            end
            for index = 1:numel(keys)
                obj.Params.(keys(index)) = values(index);
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
            if options.Step ~= 1
                error("radia:optuna:IntegerDistribution", ...
                    "Integer step values other than one are not supported yet.");
            end
            if options.Log
                if low <= 0
                    error("radia:optuna:LogBounds", ...
                        "Log integer bounds must be positive.");
                end
                key = matlab.lang.makeValidName(name);
                if isfield(obj.Params, key)
                    value = obj.Params.(key);
                    return
                end
                spec = radia.optuna.internal.DistributionCodec.integer( ...
                    low, high, true, 1);
                [isRelative, value] = obj.relativeValue(name, spec);
                if ~isRelative
                    value = obj.Study.sampleFloat(obj, name, low, high, ...
                        struct("Log", true, "Step", 1));
                    value = min(max(round(value), low), high);
                end
                obj.Params.(key) = value;
                obj.Distributions.(key) = struct( ...
                    "name", "IntDistribution", "low", low, "high", high, ...
                    "log", true, "step", 1);
                obj.Study.recordParameter(obj, name, "integer", value, ...
                    radia.optuna.internal.DistributionCodec.encode(spec));
            else
                value = obj.suggestInteger(name, low, high);
            end
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
            if ~isfinite(value) || ~isfinite(step) || step ~= floor(step)
                error("radia:optuna:Report", "value must be finite and step must be an integer.");
            end
            row = obj.IntermediateValues.Step == step;
            obj.IntermediateValues(row,:) = [];
            obj.IntermediateValues(end+1,:) = {step, value, ...
                datetime("now", "TimeZone", "local")};
            obj.IntermediateValues = sortrows(obj.IntermediateValues, "Step");
            obj.LastStep = step;
            obj.Study.recordIntermediate(obj, value, step);
        end

        function decision = shouldPrune(obj)
            decision = obj.Study.Pruner.shouldPrune(obj.Study, obj);
        end

        function decision = should_prune(obj)
            decision = obj.shouldPrune();
        end

        function prune(obj)
            obj.ensureRunning();
            obj.Study.finishTrial(obj, "PRUNED", NaN, "");
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
        function [available, value] = relativeValue(obj, name, distribution)
            key = matlab.lang.makeValidName(name);
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

        function ensureRunning(obj)
            if obj.State ~= "RUNNING"
                error("radia:optuna:TrialState", ...
                    "Trial %d is already %s.", obj.Number, obj.State);
            end
        end
    end
end
