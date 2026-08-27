classdef DistributionCodec
    %DISTRIBUTIONCODEC Persistent search-distribution metadata.

    methods (Static)
        function spec = float(low, high, logScale, step)
            spec = radia.optuna.internal.DistributionCodec.numericSpec( ...
                "float", low, high, logScale, step);
        end

        function spec = integer(low, high, logScale, step)
            spec = radia.optuna.internal.DistributionCodec.numericSpec( ...
                "integer", low, high, logScale, step);
        end

        function spec = categorical(choices)
            if isempty(choices)
                error("radia:optuna:Choices", ...
                    "Categorical choices must not be empty.");
            end
            if ischar(choices)
                choices = {choices};
            end
            spec = struct( ...
                "name", "CategoricalDistribution", ...
                "kind", "categorical", ...
                "low", NaN, ...
                "high", NaN, ...
                "log", false, ...
                "step", NaN, ...
                "choices", {choices});
        end

        function encoded = encode(spec)
            payload = struct( ...
                "schema", "radia.optuna.distribution.v1", ...
                "kind", spec.kind);
            switch spec.kind
                case {"float", "integer"}
                    payload.low = spec.low;
                    payload.high = spec.high;
                    payload.log = spec.log;
                    if isfinite(spec.step)
                        payload.step = spec.step;
                    end
                case "categorical"
                    payload.choices = spec.choices;
                otherwise
                    error("radia:optuna:DistributionKind", ...
                        "Unsupported distribution kind '%s'.", spec.kind);
            end
            encoded = string(jsonencode(payload));
        end

        function spec = decode(kind, encoded)
            kind = string(kind);
            text = string(encoded);
            decoded = [];
            decodedJSON = true;
            try
                decoded = jsondecode(char(text));
            catch
                decodedJSON = false;
            end

            if decodedJSON && isstruct(decoded) && ...
                    isfield(decoded, "schema")
                if string(decoded.schema) ~= "radia.optuna.distribution.v1"
                    error("radia:optuna:DistributionSchema", ...
                        "Unsupported distribution schema '%s'.", ...
                        string(decoded.schema));
                end
                if string(decoded.kind) ~= kind
                    error("radia:optuna:DistributionKind", ...
                        "Stored distribution kind does not match ParamTable.Kind.");
                end
                switch kind
                    case {"float", "integer"}
                        step = NaN;
                        if isfield(decoded, "step")
                            step = double(decoded.step);
                        elseif kind == "integer"
                            step = 1;
                        end
                        spec = ...
                            radia.optuna.internal.DistributionCodec.numericSpec( ...
                            kind, double(decoded.low), double(decoded.high), ...
                            logical(decoded.log), step);
                    case "categorical"
                        spec = ...
                            radia.optuna.internal.DistributionCodec.categorical( ...
                            decoded.choices);
                    otherwise
                        error("radia:optuna:DistributionKind", ...
                            "Unsupported distribution kind '%s'.", kind);
                end
                return
            end

            switch kind
                case {"float", "integer"}
                    if decodedJSON && isnumeric(decoded) && numel(decoded) == 2
                        bounds = reshape(double(decoded), 1, []);
                    else
                        bounds = sscanf(char(text), "[%f,%f]", 2).';
                    end
                    if numel(bounds) ~= 2
                        error("radia:optuna:DistributionMetadata", ...
                            "Cannot decode legacy distribution metadata '%s'.", text);
                    end
                    logScale = contains(lower(text), "log");
                    step = NaN;
                    if kind == "integer"
                        step = 1;
                    end
                    spec = ...
                        radia.optuna.internal.DistributionCodec.numericSpec( ...
                        kind, bounds(1), bounds(2), logScale, step);
                case "categorical"
                    if ~decodedJSON
                        error("radia:optuna:DistributionMetadata", ...
                            "Cannot decode categorical distribution metadata.");
                    end
                    spec = ...
                        radia.optuna.internal.DistributionCodec.categorical(decoded);
                otherwise
                    error("radia:optuna:DistributionKind", ...
                        "Unsupported distribution kind '%s'.", kind);
            end
        end

        function result = equivalent(left, right)
            if left.kind ~= right.kind
                result = false;
                return
            end
            switch left.kind
                case {"float", "integer"}
                    result = left.low == right.low && ...
                        left.high == right.high && ...
                        left.log == right.log && ...
                        isequaln(left.step, right.step);
                case "categorical"
                    result = isequal( ...
                        radia.optuna.internal.DistributionCodec.choiceTokens( ...
                        left.choices), ...
                        radia.optuna.internal.DistributionCodec.choiceTokens( ...
                        right.choices));
                otherwise
                    result = false;
            end
        end

        function result = isSingle(spec)
            if spec.kind == "categorical"
                result = isscalar(spec.choices);
            else
                result = spec.low == spec.high;
            end
        end

        function result=isSpec(value)
            required=["kind","low","high","log","step","choices"];
            result=(isa(value,"radia.optuna.BaseDistribution") && ...
                isscalar(value)) || (isstruct(value) && isscalar(value) && ...
                all(isfield(value,required)));
        end

        function spec=normalize(value)
            if isa(value,"radia.optuna.BaseDistribution") && isscalar(value)
                spec=value.toStruct();
            elseif radia.optuna.internal.DistributionCodec.isSpec(value)
                spec=value;
            else
                error("radia:optuna:Distribution", ...
                    "Value is not a Radia Optuna distribution.");
            end
        end

        function tokens = choiceTokens(choices)
            tokens = strings(numel(choices), 1);
            for index = 1:numel(choices)
                tokens(index) = ...
                    radia.optuna.internal.DistributionCodec.choiceToken( ...
                    radia.optuna.internal.DistributionCodec.choiceAt( ...
                    choices, index));
            end
        end

        function token = choiceToken(value)
            if isnumeric(value) && isscalar(value)
                token = append("numeric:", jsonencode(double(value)));
            else
                token = string(jsonencode(value));
            end
        end

        function value = choiceAt(choices, index)
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
        end
    end

    methods (Static, Access=private)
        function spec = numericSpec(kind, low, high, logScale, step)
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", ...
                    "Distribution bounds must satisfy low <= high.");
            end
            if logScale && low <= 0
                error("radia:optuna:LogBounds", ...
                    "Log-distribution bounds must be positive.");
            end
            if kind == "float" && logScale && isfinite(step)
                error("radia:optuna:Step", ...
                    "Float log distributions cannot use a finite step.");
            end
            if isfinite(step) && step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            if kind=="float"
                distributionName="FloatDistribution";
            else
                distributionName="IntDistribution";
            end
            effectiveHigh=low;
            if low~=high
                if isfinite(step)
                    effectiveHigh=low+floor((high-low)/step+1e-12)*step;
                else
                    effectiveHigh=high;
                end
            end
            spec = struct( ...
                "name", distributionName, ...
                "kind", string(kind), ...
                "low", double(low), ...
                "high", double(effectiveHigh), ...
                "log", logical(logScale), ...
                "step", double(step), ...
                "choices", {cell(0, 1)});
        end
    end
end
