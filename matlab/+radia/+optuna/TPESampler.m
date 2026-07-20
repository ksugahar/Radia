classdef TPESampler < handle
    %TPESAMPLER Tree-structured Parzen estimator sampler.
    %   This follows Optuna's user-facing sampler contract: early trials are
    %   random, then completed MATLAB table rows are split into good/bad
    %   observations and a Parzen-style density ratio chooses the next point.

    properties (SetAccess=private)
        Stream
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.25
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
    end

    methods
        function obj = TPESampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.NStartupTrials (1,1) double = 10
                options.Gamma (1,1) double = 0.25
                options.NumberOfEIChoices (1,1) double = 24
                options.PriorWeight (1,1) double = 1
            end
            if options.NStartupTrials < 0 || options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:TPEStartup", "NStartupTrials must be a nonnegative integer.");
            end
            if ~(options.Gamma > 0 && options.Gamma < 1)
                error("radia:optuna:TPEGamma", "Gamma must be between zero and one.");
            end
            obj.Stream = RandStream("mt19937ar", "Seed", double(options.Seed));
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.NumberOfEIChoices = max(1, round(options.NumberOfEIChoices));
            obj.PriorWeight = max(eps, options.PriorWeight);
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options) %#ok<INUSD>
            obj.validateBounds(low, high, options.Log, options.Step);
            [x, y] = obj.numericObservations(study, name);
            if numel(y) < obj.NStartupTrials || numel(y) < 2
                value = obj.uniform(low, high, options.Log);
                value = obj.quantize(value, low, high, options.Step);
                return;
            end

            lo = low;
            hi = high;
            if options.Log
                lo = log(low);
                hi = log(high);
            end
            if options.Log
                valid = isfinite(x) & x > 0 & x >= low & x <= high;
                x = log(x(valid));
                y = y(valid);
            else
                valid = isfinite(x) & x >= low & x <= high;
                x = x(valid);
                y = y(valid);
            end
            if numel(y) < 2
                value = obj.uniform(low, high, options.Log);
                value = obj.quantize(value, low, high, options.Step);
                return;
            end
            if study.Directions(1) == "minimize"
                [~, order] = sort(y, "ascend");
            else
                [~, order] = sort(y, "descend");
            end
            x = x(order);
            nGood = max(1, min(numel(x)-1, ceil(obj.Gamma * numel(x))));
            good = x(1:nGood);
            bad = x(nGood+1:end);
            scale = max(hi - lo, eps);
            bwGood = max(obj.bandwidth(good), scale / 1000);
            bwBad = max(obj.bandwidth(bad), scale / 1000);
            candidates = [good(:); lo + (hi-lo) * rand(obj.Stream, obj.NumberOfEIChoices, 1)];
            scores = zeros(size(candidates));
            for k = 1:numel(candidates)
                scores(k) = obj.logKde(candidates(k), good, bwGood) - ...
                    obj.logKde(candidates(k), bad, bwBad);
            end
            [~, best] = max(scores);
            if options.Log
                value = exp(candidates(best));
            else
                value = candidates(best);
            end
            value = obj.quantize(value, low, high, options.Step);
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            if low ~= floor(low) || high ~= floor(high) || low > high
                error("radia:optuna:Bounds", "Integer bounds must be finite integers with low <= high.");
            end
            value = obj.sampleFloat(study, trial, name, low, high, ...
                struct("Log", false, "Step", 1));
            value = min(max(round(value), low), high);
        end

        function value = sampleCategorical(obj, study, trial, name, choices) %#ok<INUSD>
            if isempty(choices)
                error("radia:optuna:Choices", "Categorical choices must not be empty.");
            end
            [tokens, y] = obj.categoricalObservations(study, name);
            count = numel(choices);
            if numel(y) < obj.NStartupTrials || isempty(y)
                index = 1 + floor(rand(obj.Stream, 1, 1) * count);
                value = obj.choiceAt(choices, index);
                return;
            end
            if study.Directions(1) == "minimize"
                [~, order] = sort(y, "ascend");
            else
                [~, order] = sort(y, "descend");
            end
            nGood = max(1, min(numel(y)-1, ceil(obj.Gamma * numel(y))));
            good = tokens(order(1:nGood));
            bad = tokens(order(nGood+1:end));
            scores = zeros(1, count);
            for k = 1:count
                token = obj.token( obj.choiceAt(choices, k) );
                pg = (sum(good == token) + obj.PriorWeight) / ...
                    (numel(good) + obj.PriorWeight * count);
                pb = (sum(bad == token) + obj.PriorWeight) / ...
                    (numel(bad) + obj.PriorWeight * count);
                scores(k) = log(pg) - log(pb);
            end
            [~, index] = max(scores);
            value = obj.choiceAt(choices, index);
        end

        function beforeTrial(obj, study, trial) %#ok<INUSD>
        end

        function afterTrial(obj, study, trial) %#ok<INUSD>
        end
    end

    methods (Access=private)
        function [x, y] = numericObservations(obj, study, name) %#ok<INUSL>
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & isfinite(p.ValueNumeric);
            indices = find(rows);
            x = zeros(0, 1);
            y = zeros(0, 1);
            for k = 1:numel(indices)
                row = indices(k);
                trialRow = t.TrialNumber == p.TrialNumber(row) & t.State == "COMPLETE" & ...
                    isfinite(t.Value);
                if any(trialRow)
                    x(end+1,1) = p.ValueNumeric(row); %#ok<AGROW>
                    y(end+1,1) = t.Value(find(trialRow, 1)); %#ok<AGROW>
                end
            end
        end

        function [tokens, y] = categoricalObservations(obj, study, name) %#ok<INUSL>
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & p.Kind == "categorical";
            indices = find(rows);
            tokens = strings(0,1);
            y = zeros(0,1);
            for k = 1:numel(indices)
                row = indices(k);
                trialRow = t.TrialNumber == p.TrialNumber(row) & t.State == "COMPLETE" & ...
                    isfinite(t.Value);
                if any(trialRow)
                    if isfinite(p.ValueNumeric(row))
                        tokens(end+1,1) = obj.token(p.ValueNumeric(row)); %#ok<AGROW>
                    else
                        tokens(end+1,1) = p.ValueText(row); %#ok<AGROW>
                    end
                    y(end+1,1) = t.Value(find(trialRow, 1)); %#ok<AGROW>
                end
            end
        end

        function value = uniform(obj, low, high, logScale)
            u = rand(obj.Stream, 1, 1);
            if logScale
                value = exp(log(low) + u * (log(high) - log(low)));
            else
                value = low + u * (high - low);
            end
        end

        function value = quantize(~, value, low, high, step)
            if isfinite(step)
                value = low + round((value-low)/step) * step;
            end
            value = min(max(value, low), high);
        end

        function width = bandwidth(~, values)
            if numel(values) < 2
                width = 1;
            else
                width = max(std(values), (max(values)-min(values))/max(1, sqrt(numel(values))));
            end
        end

        function value = logKde(~, point, values, width)
            z = (point - values(:)) / width;
            value = log(mean(exp(-0.5*z.^2))) - log(width) - 0.5*log(2*pi);
        end

        function validateBounds(~, low, high, logScale, step)
            if ~(isfinite(low) && isfinite(high) && low < high)
                error("radia:optuna:Bounds", "Float bounds must satisfy low < high.");
            end
            if logScale && low <= 0
                error("radia:optuna:LogBounds", "Log-uniform bounds must be positive.");
            end
            if isfinite(step) && step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
        end

        function token = token(~, value)
            if isnumeric(value) && isscalar(value)
                token = "numeric:" + string(value, "%.17g");
            else
                token = string(jsonencode(value));
            end
        end

        function value = choiceAt(~, choices, index)
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
        end
    end
end
