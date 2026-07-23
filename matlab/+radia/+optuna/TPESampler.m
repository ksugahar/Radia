classdef TPESampler < handle
    %TPESAMPLER Optuna-style tree-structured Parzen estimator sampler.
    %   The default settings and univariate Parzen construction follow the
    %   Optuna 4.9 TPESampler: ten random startup trials, a ten-percent good
    %   set capped at 25 trials, 24 expected-improvement candidates,
    %   observation-specific truncated-normal kernels, history weights,
    %   an explicit prior, and magic clipping.

    properties (SetAccess=private)
        Stream
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.1
        MaxGoodTrials (1,1) double = 25
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
        ConsiderMagicClip (1,1) logical = true
        ConsiderEndpoints (1,1) logical = false
        ConstantLiar (1,1) logical = true
        ConstraintsFcn = []
    end

    methods
        function obj = TPESampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.NStartupTrials (1,1) double = 10
                options.Gamma (1,1) double = 0.1
                options.MaxGoodTrials (1,1) double = 25
                options.NumberOfEIChoices (1,1) double = 24
                options.PriorWeight (1,1) double = 1
                options.ConsiderMagicClip (1,1) logical = true
                options.ConsiderEndpoints (1,1) logical = false
                options.ConstantLiar (1,1) logical = true
                options.ConstraintsFcn = []
            end
            if options.NStartupTrials < 0 || ...
                    options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:TPEStartup", ...
                    "NStartupTrials must be a nonnegative integer.");
            end
            if ~(options.Gamma > 0 && options.Gamma <= 1)
                error("radia:optuna:TPEGamma", ...
                    "Gamma must be greater than zero and at most one.");
            end
            if options.MaxGoodTrials < 1 || ...
                    options.MaxGoodTrials ~= floor(options.MaxGoodTrials)
                error("radia:optuna:TPEMaxGood", ...
                    "MaxGoodTrials must be a positive integer.");
            end
            if options.NumberOfEIChoices < 1 || ...
                    options.NumberOfEIChoices ~= floor(options.NumberOfEIChoices)
                error("radia:optuna:TPECandidates", ...
                    "NumberOfEIChoices must be a positive integer.");
            end
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            obj.Stream = RandStream("mt19937ar", "Seed", double(options.Seed));
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.MaxGoodTrials = options.MaxGoodTrials;
            obj.NumberOfEIChoices = options.NumberOfEIChoices;
            obj.PriorWeight = options.PriorWeight;
            obj.ConsiderMagicClip = options.ConsiderMagicClip;
            obj.ConsiderEndpoints = options.ConsiderEndpoints;
            obj.ConstantLiar = options.ConstantLiar;
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn, "function_handle")
                error("radia:optuna:ConstraintsFcn", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.ConstraintsFcn = options.ConstraintsFcn;
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options) %#ok<INUSD>
            obj.validateBounds(low, high, options.Log, options.Step);
            if low == high
                value = low;
                return
            end
            [x, y] = obj.numericObservations(study, name);
            valid = isfinite(x) & isfinite(y) & x >= low & x <= high;
            if options.Log
                valid = valid & x > 0;
            end
            x = x(valid);
            y = y(valid);
            if isempty(y) || numel(y) < obj.NStartupTrials
                value = obj.uniform(low, high, options.Log);
                value = obj.quantize(value, low, high, options.Step);
                return
            end

            [good, bad] = obj.splitObservations(x, y, study.Directions(1));
            estimatorOptions = { ...
                "Log", options.Log, ...
                "Step", options.Step, ...
                "PriorWeight", obj.PriorWeight, ...
                "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                "ConsiderEndpoints", obj.ConsiderEndpoints};
            below = radia.optuna.internal.ParzenEstimator.numerical( ...
                good, low, high, estimatorOptions{:});
            above = radia.optuna.internal.ParzenEstimator.numerical( ...
                bad, low, high, estimatorOptions{:});
            candidates = ...
                radia.optuna.internal.ParzenEstimator.sampleNumerical( ...
                below, obj.Stream, obj.NumberOfEIChoices);
            acquisition = ...
                radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                below, candidates) - ...
                radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                above, candidates);
            [~, best] = max(acquisition);
            value = candidates(best);
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            if low ~= floor(low) || high ~= floor(high) || low > high
                error("radia:optuna:Bounds", ...
                    "Integer bounds must be finite integers with low <= high.");
            end
            if low == high
                value = low;
                return
            end
            value = obj.sampleFloat(study, trial, name, low, high, ...
                struct("Log", false, "Step", 1));
            value = min(max(round(value), low), high);
        end

        function value = sampleCategorical(obj, study, trial, name, choices) %#ok<INUSD>
            if isempty(choices)
                error("radia:optuna:Choices", ...
                    "Categorical choices must not be empty.");
            end
            [tokens, y] = obj.categoricalObservations(study, name);
            choiceTokens = obj.choiceTokens(choices);
            observed = zeros(numel(tokens), 1);
            valid = isfinite(y);
            for index = 1:numel(tokens)
                match = find(choiceTokens == tokens(index), 1);
                if isempty(match)
                    valid(index) = false;
                else
                    observed(index) = match;
                end
            end
            observed = observed(valid);
            y = y(valid);
            count = numel(choiceTokens);
            if numel(y) < obj.NStartupTrials
                index = 1 + floor(rand(obj.Stream, 1, 1) * count);
                value = obj.choiceAt(choices, index);
                return
            end

            [good, bad] = obj.splitObservations( ...
                observed, y, study.Directions(1));
            below = radia.optuna.internal.ParzenEstimator.categorical( ...
                good, count, PriorWeight=obj.PriorWeight);
            above = radia.optuna.internal.ParzenEstimator.categorical( ...
                bad, count, PriorWeight=obj.PriorWeight);
            candidates = ...
                radia.optuna.internal.ParzenEstimator.sampleCategorical( ...
                below, obj.Stream, obj.NumberOfEIChoices);
            acquisition = ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                below, candidates) - ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                above, candidates);
            [~, best] = max(acquisition);
            value = obj.choiceAt(choices, candidates(best));
        end

        function values = sampleJoint(obj, study, ~, names, lows, highs, options)
            %SAMPLEJOINT Multivariate TPE with a shared mixture component.
            [observations, objectives] = obj.jointObservations(study, names, ...
                lows, highs, options.Log);
            if size(observations, 1) < obj.NStartupTrials || isempty(objectives)
                values = zeros(1, numel(names));
                for index = 1:numel(names)
                    values(index) = obj.uniform(lows(index), highs(index), ...
                        options.Log(index));
                end
                return;
            end
            [good, bad] = obj.splitJointObservations(observations, objectives, ...
                study.Directions(1));
            below = cell(1, numel(names));
            above = cell(1, numel(names));
            for index = 1:numel(names)
                estimatorOptions = {"Log", options.Log(index), "Step", NaN, ...
                    "PriorWeight", obj.PriorWeight, ...
                    "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                    "ConsiderEndpoints", obj.ConsiderEndpoints};
                below{index} = radia.optuna.internal.ParzenEstimator.numerical( ...
                    good(:,index), lows(index), highs(index), estimatorOptions{:});
                above{index} = radia.optuna.internal.ParzenEstimator.numerical( ...
                    bad(:,index), lows(index), highs(index), estimatorOptions{:});
            end
            count = obj.NumberOfEIChoices;
            components = obj.sampleComponents(below{1}.weights, count);
            candidates = zeros(count, numel(names));
            for index = 1:numel(names)
                candidates(:,index) = ...
                    radia.optuna.internal.ParzenEstimator.sampleNumericalComponents( ...
                    below{index}, obj.Stream, components);
            end
            acquisition = zeros(count, 1);
            for index = 1:numel(names)
                acquisition = acquisition + ...
                    radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                    below{index}, candidates(:,index)) - ...
                    radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                    above{index}, candidates(:,index));
            end
            [~, best] = max(acquisition);
            values = candidates(best,:);
        end

        function beforeTrial(obj, study, trial) %#ok<INUSD>
        end

        function afterTrial(obj, study, trial)
            if trial.State == "COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial, obj.ConstraintsFcn(trial));
            end
        end
    end

    methods (Access=private)
        function [good, bad] = splitObservations(obj, values, objectives, direction)
            if direction == "minimize"
                [~, order] = sort(objectives, "ascend");
            else
                [~, order] = sort(objectives, "descend");
            end
            % Match Optuna's _split_complete_trials: n_below may equal the
            % number of observations.  In that case the above density is
            % represented by its prior component only.
            nGood = min(obj.MaxGoodTrials, ceil(obj.Gamma * numel(values)));
            nGood = max(1, min(numel(values), nGood));
            isGood = false(numel(values), 1);
            isGood(order(1:nGood)) = true;
            % Keep chronological order so Optuna's history weights attach to
            % the same observations after the objective-based split.
            good = values(isGood);
            bad = values(~isGood);
        end

        function [x, y] = numericObservations(~, study, name)
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & isfinite(p.ValueNumeric);
            indices = find(rows);
            x = zeros(0, 1);
            y = zeros(0, 1);
            for index = reshape(indices, 1, [])
                trialRow = t.TrialNumber == p.TrialNumber(index) & ...
                    t.State == "COMPLETE" & isfinite(t.Value);
                if any(trialRow)
                    x(end+1, 1) = p.ValueNumeric(index); %#ok<AGROW>
                    y(end+1, 1) = t.Value(find(trialRow, 1)); %#ok<AGROW>
                end
            end
        end

        function [tokens, y] = categoricalObservations(obj, study, name)
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & p.Kind == "categorical";
            indices = find(rows);
            tokens = strings(0, 1);
            y = zeros(0, 1);
            for index = reshape(indices, 1, [])
                trialRow = t.TrialNumber == p.TrialNumber(index) & ...
                    t.State == "COMPLETE" & isfinite(t.Value);
                if any(trialRow)
                    if isfinite(p.ValueNumeric(index))
                        tokens(end+1, 1) = ...
                            obj.token(p.ValueNumeric(index)); %#ok<AGROW>
                    else
                        tokens(end+1, 1) = p.ValueText(index); %#ok<AGROW>
                    end
                    y(end+1, 1) = t.Value(find(trialRow, 1)); %#ok<AGROW>
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

        function indices = sampleComponents(obj, weights, count)
            cumulative = cumsum(weights(:));
            cumulative(end) = 1;
            u = rand(obj.Stream, count, 1);
            indices = 1 + sum(cumulative.' < u, 2);
        end

        function [x, y] = jointObservations(obj, study, names, lows, highs, logs)
            states = study.TrialTable.State;
            usable = states == "COMPLETE" | states == "PRUNED" | ...
                (states == "RUNNING" & obj.ConstantLiar);
            trialNumbers = study.TrialTable.TrialNumber(usable);
            finished = study.TrialTable.Value(usable);
            finiteFinished = finished(isfinite(finished));
            if isempty(finiteFinished)
                liar = 0;
            elseif study.Directions(1) == "minimize"
                liar = max(finiteFinished);
            else
                liar = min(finiteFinished);
            end
            x = zeros(0, numel(names));
            y = zeros(0, 1);
            for number = reshape(trialNumbers, 1, [])
                row = study.TrialTable.TrialNumber == number;
                state = study.TrialTable.State(find(row,1));
                values = NaN(1, numel(names));
                valid = true;
                for index = 1:numel(names)
                    p = study.ParamTable.TrialNumber == number & ...
                        study.ParamTable.Name == names(index) & ...
                        study.ParamTable.Kind == "float";
                    if ~any(p)
                        valid = false;
                        break;
                    end
                    values(index) = study.ParamTable.ValueNumeric(find(p,1));
                    valid = valid && values(index) >= lows(index) && ...
                        values(index) <= highs(index) && ...
                        (~logs(index) || values(index) > 0);
                end
                if state == "COMPLETE"
                    objectiveValue = study.TrialTable.Value(find(row,1));
                elseif state == "PRUNED"
                    intermediate = study.IntermediateTable.TrialNumber == number;
                    if any(intermediate)
                        objectiveValue = study.IntermediateTable.Value( ...
                            find(intermediate,1,'last'));
                    else
                        valid = false;
                        objectiveValue = NaN;
                    end
                else
                    objectiveValue = liar;
                end
                if valid
                    x(end+1,:) = values; %#ok<AGROW>
                    y(end+1,1) = objectiveValue; %#ok<AGROW>
                end
            end
        end

        function [good, bad] = splitJointObservations(obj, values, objectives, direction)
            if direction == "minimize"
                [~, order] = sort(objectives, "ascend");
            else
                [~, order] = sort(objectives, "descend");
            end
            nGood = min(obj.MaxGoodTrials, ceil(obj.Gamma * size(values,1)));
            nGood = max(1, min(size(values,1), nGood));
            mask = false(size(values,1),1);
            mask(order(1:nGood)) = true;
            good = values(mask,:);
            bad = values(~mask,:);
        end

        function value = quantize(~, value, low, high, step)
            if isfinite(step)
                value = low + round((value - low) / step) * step;
            end
            value = min(max(value, low), high);
        end

        function validateBounds(~, low, high, logScale, step)
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", ...
                    "Float bounds must satisfy low <= high.");
            end
            if logScale && low <= 0
                error("radia:optuna:LogBounds", ...
                    "Log-uniform bounds must be positive.");
            end
            if isfinite(step) && step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            if logScale && isfinite(step) && low - step / 2 <= 0
                error("radia:optuna:LogBounds", ...
                    "Expanded log-distribution support must be positive.");
            end
        end

        function tokens = choiceTokens(obj, choices)
            tokens = strings(numel(choices), 1);
            for index = 1:numel(choices)
                tokens(index) = obj.token(obj.choiceAt(choices, index));
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
