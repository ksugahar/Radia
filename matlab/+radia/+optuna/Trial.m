classdef Trial < handle
    %TRIAL One define-by-run evaluation in a table-backed Study.

    properties (SetAccess=private)
        Study
        Number (1,1) double
        State (1,1) string = "RUNNING"
        Value (1,1) double = NaN
        Values double = NaN
        Params struct = struct()
        IntermediateValues table
        UserAttrs struct = struct()
        StartTime datetime
        EndTime datetime = NaT
        ErrorMessage (1,1) string = ""
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
            value = obj.Study.sampleFloat(obj, name, low, high, ...
                struct("Log", options.Log, "Step", options.Step));
            obj.Params.(key) = value;
            obj.Study.recordParameter(obj, name, "float", value, ...
                sprintf("[%g,%g]", low, high));
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
            value = obj.Study.sampleInteger(obj, name, low, high);
            obj.Params.(key) = value;
            obj.Study.recordParameter(obj, name, "integer", value, ...
                sprintf("[%g,%g]", low, high));
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
            value = obj.Study.sampleCategorical(obj, name, choices);
            obj.Params.(key) = value;
            obj.Study.recordParameter(obj, name, "categorical", value, ...
                string(jsonencode(choices)));
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
            if options.Step ~= 1 || options.Log
                error("radia:optuna:IntegerDistribution", ...
                    "This MATLAB integer distribution currently supports unit steps and linear sampling.");
            end
            value = obj.suggestInteger(name, low, high);
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
    end

    methods (Access=private)
        function ensureRunning(obj)
            if obj.State ~= "RUNNING"
                error("radia:optuna:TrialState", ...
                    "Trial %d is already %s.", obj.Number, obj.State);
            end
        end
    end
end
