classdef CmaEsSampler < handle
    %CMAESSAMPLER Stateful CMA-ES-compatible sampler for numeric variables.
    %   The public behavior follows Optuna's sampler lifecycle.  Continuous
    %   variables are optimized jointly; integer variables are rounded and
    %   categorical variables use a seeded random fallback.

    properties (SetAccess=private)
        Stream
        NStartupTrials (1,1) double = 1
        PopulationSize (1,1) double = 0
        Sigma (1,1) double = 0.3
    end

    properties (Access=private)
        DimensionNames string = strings(1,0)
        LogDimensions logical = false(1,0)
        LowerBounds double = zeros(1,0)
        UpperBounds double = zeros(1,0)
        Mean double = zeros(1,0)
        Covariance double = zeros(0,0)
        ActiveTrialNumber (1,1) double = -1
        ActiveCandidate double = zeros(1,0)
        ActiveValues struct = struct()
        HistoryX double = zeros(0,0)
        HistoryY double = zeros(0,1)
    end

    methods
        function obj = CmaEsSampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.NStartupTrials (1,1) double = 1
                options.PopulationSize (1,1) double = 0
                options.Sigma (1,1) double = 0.3
            end
            if options.NStartupTrials < 0 || options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:CMAStartup", "NStartupTrials must be a nonnegative integer.");
            end
            if ~(options.Sigma > 0 && isfinite(options.Sigma))
                error("radia:optuna:CMASigma", "Sigma must be positive and finite.");
            end
            obj.Stream = RandStream("mt19937ar", "Seed", double(options.Seed));
            obj.NStartupTrials = options.NStartupTrials;
            obj.PopulationSize = options.PopulationSize;
            obj.Sigma = options.Sigma;
        end

        function beforeTrial(obj, study, trial) %#ok<INUSD>
            obj.ActiveTrialNumber = trial.Number;
            obj.ActiveCandidate = zeros(1,0);
            obj.ActiveValues = struct();
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options) %#ok<INUSD>
            if ~(isfinite(low) && isfinite(high) && low < high)
                error("radia:optuna:Bounds", "Float bounds must satisfy low < high.");
            end
            if options.Log && low <= 0
                error("radia:optuna:LogBounds", "Log-uniform bounds must be positive.");
            end
            key = matlab.lang.makeValidName(name);
            if isfield(obj.ActiveValues, key)
                value = obj.ActiveValues.(key);
                return;
            end
            index = obj.ensureDimension(name, low, high, options.Log);
            if trial.Number < obj.NStartupTrials
                u = rand(obj.Stream, 1, 1);
                if options.Log
                    value = exp(log(low) + u * (log(high)-log(low)));
                else
                    value = low + u * (high-low);
                end
            else
                obj.ensureCandidate();
                value = obj.ActiveCandidate(index);
                if options.Log
                    value = exp(value);
                end
                value = min(max(value, low), high);
            end
            if isfinite(options.Step)
                if options.Step <= 0
                    error("radia:optuna:Step", "Step must be positive.");
                end
                value = low + round((value-low)/options.Step)*options.Step;
                value = min(max(value, low), high);
            end
            obj.ActiveValues.(key) = value;
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            if low ~= floor(low) || high ~= floor(high) || low > high
                error("radia:optuna:Bounds", "Integer bounds must be finite integers with low <= high.");
            end
            value = obj.sampleFloat(study, trial, name, low, high, ...
                struct("Log", false, "Step", 1));
            value = min(max(round(value), low), high);
            obj.ActiveValues.(matlab.lang.makeValidName(name)) = value;
        end

        function value = sampleCategorical(obj, study, trial, name, choices) %#ok<INUSD>
            if isempty(choices)
                error("radia:optuna:Choices", "Categorical choices must not be empty.");
            end
            key = matlab.lang.makeValidName(name);
            if isfield(obj.ActiveValues, key)
                value = obj.ActiveValues.(key);
                return;
            end
            count = numel(choices);
            index = 1 + floor(rand(obj.Stream, 1, 1) * count);
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
            obj.ActiveValues.(key) = value;
        end

        function afterTrial(obj, study, trial)
            if trial.State ~= "COMPLETE" || isempty(obj.DimensionNames)
                return;
            end
            row = NaN(1, numel(obj.DimensionNames));
            for k = 1:numel(obj.DimensionNames)
                key = matlab.lang.makeValidName(obj.DimensionNames(k));
                if ~isfield(trial.Params, key)
                    return;
                end
                value = trial.Params.(key);
                if obj.LogDimensions(k)
                    if ~(isfinite(value) && value > 0)
                        return;
                    end
                    row(k) = log(value);
                else
                    row(k) = value;
                end
            end
            if any(~isfinite(row))
                return;
            end
            obj.HistoryX(end+1,:) = row;
            obj.HistoryY(end+1,1) = trial.Value;
            obj.updateDistribution(study);
        end
    end

    methods (Access=private)
        function index = ensureDimension(obj, name, low, high, isLog)
            internalLow = low;
            internalHigh = high;
            if isLog
                internalLow = log(low);
                internalHigh = log(high);
            end
            index = find(obj.DimensionNames == string(name), 1);
            if ~isempty(index)
                if obj.LogDimensions(index) ~= isLog || ...
                        abs(obj.LowerBounds(index) - internalLow) > eps(max(1, abs(internalLow))) || ...
                        abs(obj.UpperBounds(index) - internalHigh) > eps(max(1, abs(internalHigh)))
                    error("radia:optuna:DistributionChanged", ...
                        "Distribution for parameter '%s' changed during the study.", name);
                end
                return;
            end
            obj.DimensionNames(end+1) = string(name);
            obj.LogDimensions(end+1) = isLog;
            obj.LowerBounds(end+1) = internalLow;
            obj.UpperBounds(end+1) = internalHigh;
            obj.Mean(end+1) = 0.5*(internalLow+internalHigh);
            d = numel(obj.Mean);
            if ~isempty(obj.HistoryX) && size(obj.HistoryX, 2) < d
                obj.HistoryX(:, end+1:d) = NaN;
            end
            if d == 1
                obj.Covariance = 1;
            else
                old = obj.Covariance;
                obj.Covariance = blkdiag(old, 1);
            end
            index = d;
        end

        function ensureCandidate(obj)
            d = numel(obj.Mean);
            if isempty(obj.ActiveCandidate)
                [L, flag] = chol(obj.Covariance + 1e-10*eye(d), "lower");
                if flag ~= 0
                    L = eye(d);
                end
                z = randn(obj.Stream, d, 1);
                widths = max(obj.UpperBounds-obj.LowerBounds, eps).';
                obj.ActiveCandidate = obj.Mean + ...
                    (obj.Sigma * (widths .* (L*z))).';
            elseif numel(obj.ActiveCandidate) < d
                newIndex = numel(obj.ActiveCandidate)+1;
                width = max(obj.UpperBounds(newIndex)-obj.LowerBounds(newIndex), eps);
                obj.ActiveCandidate(newIndex) = obj.Mean(newIndex) + ...
                    obj.Sigma*width*randn(obj.Stream, 1, 1);
            end
            for k = 1:d
                width = max(obj.UpperBounds(k)-obj.LowerBounds(k), eps);
                obj.ActiveCandidate(k) = min(max(obj.ActiveCandidate(k), ...
                    obj.LowerBounds(k)-obj.Sigma*width), ...
                    obj.UpperBounds(k)+obj.Sigma*width);
            end
        end

        function updateDistribution(obj, study)
            if isempty(obj.HistoryX)
                return;
            end
            valid = all(isfinite(obj.HistoryX), 2) & isfinite(obj.HistoryY);
            if ~any(valid)
                return;
            end
            historyX = obj.HistoryX(valid, :);
            historyY = obj.HistoryY(valid);
            if study.Directions(1) == "minimize"
                [~, order] = sort(historyY, "ascend");
            else
                [~, order] = sort(historyY, "descend");
            end
            d = size(historyX, 2);
            lambda = obj.PopulationSize;
            if lambda <= 0
                lambda = max(4, 4 + floor(3*log(max(2,d))));
            end
            mu = max(1, min(numel(order), floor(lambda/2)));
            weights = log(mu + 0.5) - log((1:mu).');
            weights = weights / sum(weights);
            selected = historyX(order(1:mu), :);
            newMean = (weights.' * selected);
            centered = selected - newMean;
            newCov = zeros(d,d);
            for k = 1:mu
                newCov = newCov + weights(k) * (centered(k,:).'*centered(k,:));
            end
            scales = max(obj.UpperBounds-obj.LowerBounds, eps);
            newCov = newCov ./ (scales.'*scales);
            newCov = 0.8*obj.Covariance + 0.2*(newCov + 1e-8*eye(d));
            obj.Mean = min(max(newMean, obj.LowerBounds), obj.UpperBounds);
            obj.Covariance = (newCov + newCov.')/2;
            spread = sqrt(max(diag(obj.Covariance), 1e-8)).';
            obj.Sigma = min(1.0, max(0.02, 0.9*obj.Sigma + 0.1*mean(spread)));
        end
    end
end
