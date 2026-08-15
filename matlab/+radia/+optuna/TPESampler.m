classdef TPESampler < handle
    %TPESAMPLER Optuna-style tree-structured Parzen estimator sampler.
    %   The default settings and univariate Parzen construction follow the
    %   Optuna 4.9 TPESampler: ten random startup trials, a ten-percent good
    %   set capped at 25 trials, 24 expected-improvement candidates,
    %   observation-specific truncated-normal kernels, history weights,
    %   an explicit prior, and magic clipping.
    %   Constraint-aware splits use c <= 0 as feasible, rank infeasible
    %   trials by total positive violation, and rank missing/nonfinite
    %   constraint records last.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.1
        MaxGoodTrials (1,1) double = 25
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
        ConsiderMagicClip (1,1) logical = true
        ConsiderEndpoints (1,1) logical = false
        Multivariate (1,1) logical = false
        ConstantLiar (1,1) logical = false
        ConstraintsFcn = []
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.tpe-sampler-state.v1"
        SamplerName = "tpe"
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
                options.Multivariate (1,1) logical = false
                options.ConstantLiar (1,1) logical = false
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
            obj.Seed = double(options.Seed);
            obj.Stream = RandStream("mt19937ar", "Seed", obj.Seed);
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.MaxGoodTrials = options.MaxGoodTrials;
            obj.NumberOfEIChoices = options.NumberOfEIChoices;
            obj.PriorWeight = options.PriorWeight;
            obj.ConsiderMagicClip = options.ConsiderMagicClip;
            obj.ConsiderEndpoints = options.ConsiderEndpoints;
            obj.Multivariate = options.Multivariate;
            obj.ConstantLiar = options.ConstantLiar;
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn, "function_handle")
                error("radia:optuna:ConstraintsFcn", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.ConstraintsFcn = options.ConstraintsFcn;
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options)
            obj.attach(study);
            obj.validateBounds(low, high, options.Log, options.Step);
            if low == high
                value = low;
                return
            end
            [x, y, trialNumbers, pending] = ...
                obj.numericObservations(study, name);
            valid = isfinite(x) & isfinite(y) & x >= low & x <= high;
            if options.Log
                valid = valid & x > 0;
            end
            x = x(valid);
            y = y(valid);
            trialNumbers = trialNumbers(valid);
            pending = pending(valid);
            finishedCount = sum(~pending);
            if finishedCount == 0 || finishedCount < obj.NStartupTrials
                value = obj.randomNumerical( ...
                    low, high, options.Log, options.Step);
                obj.recordState(study, trial.Number);
                return
            end

            [good, bad] = obj.splitObservations( ...
                x, y, study.Directions(1), study, trialNumbers, pending);
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
            obj.recordState(study, trial.Number);
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

        function value = sampleCategorical(obj, study, trial, name, choices)
            obj.attach(study);
            if isempty(choices)
                error("radia:optuna:Choices", ...
                    "Categorical choices must not be empty.");
            end
            [tokens, y, trialNumbers, pending] = ...
                obj.categoricalObservations(study, name);
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
            trialNumbers = trialNumbers(valid);
            pending = pending(valid);
            count = numel(choiceTokens);
            finishedCount = sum(~pending);
            if finishedCount == 0 || finishedCount < obj.NStartupTrials
                index = 1 + floor(rand(obj.Stream, 1, 1) * count);
                value = obj.choiceAt(choices, index);
                obj.recordState(study, trial.Number);
                return
            end

            [good, bad] = obj.splitObservations( ...
                observed, y, study.Directions(1), study, ...
                trialNumbers, pending);
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
            obj.recordState(study, trial.Number);
        end

        function values = sampleJoint(obj, study, trial, names, lows, highs, options)
            %SAMPLEJOINT Multivariate TPE with a shared mixture component.
            obj.attach(study);
            template = struct( ...
                "name", "", ...
                "distribution", ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    0, 1, false, NaN));
            searchSpace = repmat(template, 1, numel(names));
            for index = 1:numel(names)
                searchSpace(index).name = names(index);
                searchSpace(index).distribution = ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    lows(index), highs(index), options.Log(index), NaN);
            end
            sampled = obj.sampleRelativeSpace(study, searchSpace);
            values = zeros(1, numel(sampled));
            for index = 1:numel(sampled)
                values(index) = double(sampled{index});
            end
            obj.recordState(study, trial.Number);
        end

        function searchSpace = inferRelativeSearchSpace(obj, study, trial) %#ok<INUSD>
            if ~obj.Multivariate
                searchSpace = obj.emptySearchSpace();
                return
            end
            searchSpace = obj.intersectionSearchSpace(study);
        end

        function searchSpace = infer_relative_search_space(obj, study, trial)
            if nargin < 3
                trial = [];
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
        end

        function beforeTrial(obj, study, trial)
            obj.attach(study);
            if ~obj.Multivariate || numel(study.Directions) ~= 1
                return
            end
            finished = study.TrialTable.State == "COMPLETE" | ...
                study.TrialTable.State == "PRUNED";
            if sum(finished) < obj.NStartupTrials
                return
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
            if isempty(searchSpace)
                return
            end
            values = obj.sampleRelativeSpace(study, searchSpace);
            trial.setRelativeParameters(searchSpace, values);
            obj.recordState(study, trial.Number);
        end

        function afterTrial(obj, study, trial)
            if trial.State == "COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial, obj.ConstraintsFcn(trial));
            end
        end
    end

    methods (Access=private)
        function [good, bad] = splitObservations(obj, values, objectives, ...
                direction, study, trialNumbers, pending)
            % Match Optuna's _split_complete_trials: n_below may equal the
            % number of observations.  In that case the above density is
            % represented by its prior component only.
            finishedCount = sum(~pending);
            nGood = min(obj.MaxGoodTrials, ...
                ceil(obj.Gamma * finishedCount));
            nGood = max(1, min(finishedCount, nGood));
            order = obj.rankObservations( ...
                objectives, direction, study, trialNumbers, pending);
            isGood = false(numel(values), 1);
            isGood(order(1:nGood)) = true;
            % Keep chronological order so Optuna's history weights attach to
            % the same observations after the objective-based split.
            good = values(isGood);
            bad = values(~isGood);
        end

        function [x, y, trialNumbers, pending] = ...
                numericObservations(obj, study, name)
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & isfinite(p.ValueNumeric);
            indices = find(rows);
            x = zeros(0, 1);
            y = zeros(0, 1);
            trialNumbers = zeros(0, 1);
            pending = false(0, 1);
            liar = obj.liarObjective(study);
            for index = reshape(indices, 1, [])
                trialRow = find(t.TrialNumber == p.TrialNumber(index), 1);
                if isempty(trialRow)
                    continue
                end
                state = t.State(trialRow);
                if state == "COMPLETE" && isfinite(t.Value(trialRow))
                    objective = t.Value(trialRow);
                    isPending = false;
                elseif state == "RUNNING" && obj.ConstantLiar
                    objective = liar;
                    isPending = true;
                else
                    continue
                end
                x(end+1, 1) = p.ValueNumeric(index); %#ok<AGROW>
                y(end+1, 1) = objective; %#ok<AGROW>
                trialNumbers(end+1, 1) = p.TrialNumber(index); %#ok<AGROW>
                pending(end+1, 1) = isPending; %#ok<AGROW>
            end
        end

        function [tokens, y, trialNumbers, pending] = ...
                categoricalObservations(obj, study, name)
            p = study.ParamTable;
            t = study.TrialTable;
            rows = p.Name == string(name) & p.Kind == "categorical";
            indices = find(rows);
            tokens = strings(0, 1);
            y = zeros(0, 1);
            trialNumbers = zeros(0, 1);
            pending = false(0, 1);
            liar = obj.liarObjective(study);
            for index = reshape(indices, 1, [])
                trialRow = find(t.TrialNumber == p.TrialNumber(index), 1);
                if isempty(trialRow)
                    continue
                end
                state = t.State(trialRow);
                if state == "COMPLETE" && isfinite(t.Value(trialRow))
                    objective = t.Value(trialRow);
                    isPending = false;
                elseif state == "RUNNING" && obj.ConstantLiar
                    objective = liar;
                    isPending = true;
                else
                    continue
                end
                if isfinite(p.ValueNumeric(index))
                    tokens(end+1, 1) = ...
                        obj.token(p.ValueNumeric(index)); %#ok<AGROW>
                else
                    tokens(end+1, 1) = p.ValueText(index); %#ok<AGROW>
                end
                y(end+1, 1) = objective; %#ok<AGROW>
                trialNumbers(end+1, 1) = p.TrialNumber(index); %#ok<AGROW>
                pending(end+1, 1) = isPending; %#ok<AGROW>
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

        function searchSpace = intersectionSearchSpace(~, study)
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.calculate( ...
                study, IncludePruned=true);
        end

        function searchSpace = emptySearchSpace(~)
            template = struct( ...
                "name", "", ...
                "distribution", ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    0, 1, false, NaN));
            searchSpace = reshape(template([]), 0, 1);
        end

        function values = sampleRelativeSpace(obj, study, searchSpace)
            [observations, objectives, trialNumbers, pending] = ...
                obj.relativeObservations(study, searchSpace);
            finishedCount = sum(~pending);
            if finishedCount == 0 || ...
                    finishedCount < obj.NStartupTrials || ...
                    isempty(objectives)
                values = cell(1, numel(searchSpace));
                for index = 1:numel(searchSpace)
                    values{index} = obj.randomRelativeValue( ...
                        searchSpace(index).distribution);
                end
                return
            end

            [good, bad] = obj.splitJointObservations( ...
                observations, objectives, study.Directions(1), study, ...
                trialNumbers, pending);
            dimension = numel(searchSpace);
            below = cell(1, dimension);
            above = cell(1, dimension);
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    choiceCount = numel(distribution.choices);
                    below{index} = ...
                        radia.optuna.internal.ParzenEstimator.categorical( ...
                        good(:,index), choiceCount, ...
                        PriorWeight=obj.PriorWeight);
                    above{index} = ...
                        radia.optuna.internal.ParzenEstimator.categorical( ...
                        bad(:,index), choiceCount, ...
                        PriorWeight=obj.PriorWeight);
                else
                    estimatorOptions = { ...
                        "Log", distribution.log, ...
                        "Step", distribution.step, ...
                        "PriorWeight", obj.PriorWeight, ...
                        "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                        "ConsiderEndpoints", obj.ConsiderEndpoints, ...
                        "MultivariateDimension", dimension};
                    below{index} = ...
                        radia.optuna.internal.ParzenEstimator.numerical( ...
                        good(:,index), distribution.low, distribution.high, ...
                        estimatorOptions{:});
                    above{index} = ...
                        radia.optuna.internal.ParzenEstimator.numerical( ...
                        bad(:,index), distribution.low, distribution.high, ...
                        estimatorOptions{:});
                end
            end

            count = obj.NumberOfEIChoices;
            components = obj.sampleComponents(below{1}.weights, count);
            candidates = zeros(count, dimension);
            acquisition = zeros(count, 1);
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    candidates(:,index) = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        sampleCategoricalComponents( ...
                        below{index}, obj.Stream, components);
                    acquisition = acquisition + ...
                        radia.optuna.internal.ParzenEstimator. ...
                        logPdfCategorical(below{index}, candidates(:,index)) - ...
                        radia.optuna.internal.ParzenEstimator. ...
                        logPdfCategorical(above{index}, candidates(:,index));
                else
                    candidates(:,index) = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        sampleNumericalComponents( ...
                        below{index}, obj.Stream, components);
                    acquisition = acquisition + ...
                        radia.optuna.internal.ParzenEstimator. ...
                        logPdfNumerical(below{index}, candidates(:,index)) - ...
                        radia.optuna.internal.ParzenEstimator. ...
                        logPdfNumerical(above{index}, candidates(:,index));
                end
            end
            [~, best] = max(acquisition);
            values = cell(1, dimension);
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    values{index} = ...
                        radia.optuna.internal.DistributionCodec.choiceAt( ...
                        distribution.choices, candidates(best,index));
                else
                    values{index} = candidates(best,index);
                end
            end
        end

        function [x, y, observationTrialNumbers, pending] = ...
                relativeObservations(obj, study, searchSpace)
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
            x = zeros(0, numel(searchSpace));
            y = zeros(0, 1);
            observationTrialNumbers = zeros(0, 1);
            pending = false(0, 1);
            for number = reshape(trialNumbers, 1, [])
                row = study.TrialTable.TrialNumber == number;
                state = study.TrialTable.State(find(row,1));
                values = NaN(1, numel(searchSpace));
                valid = true;
                for index = 1:numel(searchSpace)
                    p = study.ParamTable.TrialNumber == number & ...
                        study.ParamTable.Name == searchSpace(index).name;
                    if ~any(p)
                        valid = false;
                        break;
                    end
                    parameterRow = find(p, 1);
                    stored = ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        study.ParamTable.Kind(parameterRow), ...
                        study.ParamTable.Distribution(parameterRow));
                    valid = valid && ...
                        radia.optuna.internal.DistributionCodec.equivalent( ...
                        searchSpace(index).distribution, stored);
                    if ~valid
                        break
                    end
                    [validValue, values(index)] = obj.parameterInternalValue( ...
                        study.ParamTable(parameterRow,:), ...
                        searchSpace(index).distribution);
                    valid = valid && validValue;
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
                    observationTrialNumbers(end+1,1) = number; %#ok<AGROW>
                    pending(end+1,1) = state == "RUNNING"; %#ok<AGROW>
                end
            end
        end

        function [valid, value] = parameterInternalValue(obj, row, distribution)
            if distribution.kind == "categorical"
                if isfinite(row.ValueNumeric)
                    token = obj.token(row.ValueNumeric);
                else
                    token = row.ValueText;
                end
                tokens = ...
                    radia.optuna.internal.DistributionCodec.choiceTokens( ...
                    distribution.choices);
                match = find(tokens == token, 1);
                valid = ~isempty(match);
                if valid
                    value = match;
                else
                    value = NaN;
                end
                return
            end
            value = row.ValueNumeric;
            valid = isfinite(value) && value >= distribution.low && ...
                value <= distribution.high && ...
                (~distribution.log || value > 0);
        end

        function value = randomRelativeValue(obj, distribution)
            if distribution.kind == "categorical"
                index = 1 + floor(rand(obj.Stream, 1, 1) * ...
                    numel(distribution.choices));
                value = ...
                    radia.optuna.internal.DistributionCodec.choiceAt( ...
                    distribution.choices, index);
                return
            end
            value = obj.uniform( ...
                distribution.low, distribution.high, distribution.log);
            value = obj.quantize(value, distribution.low, ...
                distribution.high, distribution.step);
            if distribution.kind == "integer"
                value = round(value);
            end
        end

        function [good, bad] = splitJointObservations(obj, values, ...
                objectives, direction, study, trialNumbers, pending)
            finishedCount = sum(~pending);
            nGood = min(obj.MaxGoodTrials, ...
                ceil(obj.Gamma * finishedCount));
            nGood = max(1, min(finishedCount, nGood));
            order = obj.rankObservations( ...
                objectives, direction, study, trialNumbers, pending);
            mask = false(size(values,1),1);
            mask(order(1:nGood)) = true;
            good = values(mask,:);
            bad = values(~mask,:);
        end

        function order = rankObservations(obj, objectives, direction, ...
                study, trialNumbers, pending)
            % Feasible completed/pruned trials rank before infeasible ones;
            % infeasible trials rank by total positive violation. Missing or
            % nonfinite constraint data is deliberately worst (Inf). Pending
            % constant-liar observations can shape only the above density.
            violations = obj.constraintViolations(study, trialNumbers);
            objectiveRank = reshape(double(objectives), [], 1);
            if direction == "maximize"
                objectiveRank = -objectiveRank;
            end
            sequence = (1:numel(objectiveRank))';
            rankKeys = [double(reshape(pending, [], 1)), ...
                double(violations > 0), violations, objectiveRank, sequence];
            [~, order] = sortrows(rankKeys, 1:size(rankKeys, 2));
        end

        function violations = constraintViolations(obj, study, trialNumbers)
            trialNumbers = reshape(double(trialNumbers), [], 1);
            violations = zeros(size(trialNumbers));
            constraintsEnabled = ~isempty(obj.ConstraintsFcn) || ...
                study.hasConstraintRecords();
            if ~constraintsEnabled
                return
            end
            for index = 1:numel(trialNumbers)
                [present,constraintValues] = ...
                    study.constraintRecord(trialNumbers(index));
                if ~present
                    violations(index) = Inf;
                    continue
                end
                violations(index) = sum(max(constraintValues, 0));
            end
        end

        function value = liarObjective(~, study)
            complete = study.TrialTable.State == "COMPLETE" & ...
                isfinite(study.TrialTable.Value);
            finished = study.TrialTable.Value(complete);
            if isempty(finished)
                value = 0;
            elseif study.Directions(1) == "minimize"
                value = max(finished);
            else
                value = min(finished);
            end
        end

        function value = randomNumerical(obj, low, high, logScale, step)
            if isfinite(step) && ~logScale
                count = floor((high - low) / step + 1e-12) + 1;
                index = floor(rand(obj.Stream, 1, 1) * count);
                value = low + index * step;
                value = min(value, high);
                return
            end
            value = obj.uniform(low, high, logScale);
            value = obj.quantize(value, low, high, step);
        end

        function attach(obj, study)
            changed = isempty(obj.AttachedStudy) || ...
                ~isequal(obj.AttachedStudy, study);
            if changed
                obj.AttachedStudy = study;
                obj.Stream = RandStream("mt19937ar", "Seed", obj.Seed);
                obj.Restored = false;
            end
            if obj.Restored
                return
            end
            state = study.samplerState(obj.SamplerName, obj.StateSchema);
            if ~isempty(state)
                obj.restoreState(state);
            end
            obj.Restored = true;
        end

        function restoreState(obj, state)
            required = ["schema", "seed", "random_state"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state, required)) || ...
                    string(state.schema) ~= obj.StateSchema
                error("radia:optuna:TPEState", ...
                    "Stored TPE sampler state is invalid or unsupported.");
            end
            if ~isnumeric(state.seed) || ~isscalar(state.seed) || ...
                    ~isfinite(double(state.seed)) || ...
                    double(state.seed) ~= obj.Seed
                error("radia:optuna:TPEStateSeed", ...
                    "Stored TPE sampler seed (%g) does not match the " + ...
                    "configured seed (%g).", double(state.seed), obj.Seed);
            end
            try
                obj.Stream.State = state.random_state;
            catch exception
                error("radia:optuna:TPEState", ...
                    "Stored TPE random state is invalid: %s", ...
                    exception.message);
            end
        end

        function recordState(obj, study, trialNumber)
            state = struct( ...
                "schema", obj.StateSchema, ...
                "seed", obj.Seed, ...
                "random_state", obj.Stream.State);
            generation = sum(study.TrialTable.State == "COMPLETE");
            study.recordSamplerState(obj.SamplerName, obj.StateSchema, ...
                trialNumber, generation, state);
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
