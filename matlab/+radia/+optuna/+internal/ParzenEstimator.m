classdef ParzenEstimator
    %PARZENESTIMATOR Optuna-style mixture estimator for MATLAB TPE.
    %   Numerical kernels use observation-specific truncated-normal widths,
    %   an explicit domain prior, magic clipping, and Optuna's history
    %   weights. Categorical kernels use one smoothed distribution per
    %   observation plus a uniform prior component.

    methods (Static)
        function estimator = numerical(observations, low, high, options)
            arguments
                observations double
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = NaN
                options.PriorWeight (1,1) double = 1
                options.ConsiderMagicClip (1,1) logical = true
                options.ConsiderEndpoints (1,1) logical = false
            end
            observations = reshape(double(observations), [], 1);
            if any(~isfinite(observations))
                error("radia:optuna:ParzenObservations", ...
                    "Parzen observations must be finite.");
            end
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            if ~(isfinite(low) && isfinite(high) && low < high)
                error("radia:optuna:Bounds", ...
                    "Float bounds must satisfy low < high.");
            end
            isDiscrete = isfinite(options.Step);
            if isDiscrete && options.Step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end

            supportLow = low;
            supportHigh = high;
            if isDiscrete
                supportLow = low - options.Step / 2;
                supportHigh = high + options.Step / 2;
            end
            if options.Log
                if supportLow <= 0
                    error("radia:optuna:LogBounds", ...
                        "Log-distribution support must be positive.");
                end
                mus = log(observations);
                internalLow = log(supportLow);
                internalHigh = log(supportHigh);
            else
                mus = observations;
                internalLow = supportLow;
                internalHigh = supportHigh;
            end

            n = numel(mus);
            sigmas = radia.optuna.internal.ParzenEstimator.sigmas( ...
                mus, internalLow, internalHigh, ...
                options.ConsiderMagicClip, options.ConsiderEndpoints);
            mus = [mus; 0.5 * (internalLow + internalHigh)];
            sigmas = [sigmas; internalHigh - internalLow];
            weights = radia.optuna.internal.ParzenEstimator.mixtureWeights( ...
                n, options.PriorWeight);

            estimator = struct( ...
                "kind", "numerical", ...
                "weights", weights, ...
                "mu", mus, ...
                "sigma", sigmas, ...
                "low", low, ...
                "high", high, ...
                "support_low", supportLow, ...
                "support_high", supportHigh, ...
                "internal_low", internalLow, ...
                "internal_high", internalHigh, ...
                "log", options.Log, ...
                "step", options.Step);
        end

        function estimator = categorical(observedIndices, nChoices, options)
            arguments
                observedIndices double
                nChoices (1,1) double {mustBeInteger, mustBePositive}
                options.PriorWeight (1,1) double = 1
            end
            observedIndices = reshape(double(observedIndices), [], 1);
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            if any(observedIndices ~= floor(observedIndices)) || ...
                    any(observedIndices < 1 | observedIndices > nChoices)
                error("radia:optuna:ParzenCategory", ...
                    "Categorical observations must index the available choices.");
            end

            n = numel(observedIndices);
            if n == 0
                probabilities = ones(1, nChoices) / nChoices;
            else
                nKernels = n + 1;
                probabilities = repmat(options.PriorWeight / nKernels, ...
                    nKernels, nChoices);
                for index = 1:n
                    probabilities(index, observedIndices(index)) = ...
                        probabilities(index, observedIndices(index)) + 1;
                end
                rowSums = sum(probabilities, 2);
                zeroRows = rowSums == 0;
                probabilities(~zeroRows, :) = ...
                    probabilities(~zeroRows, :) ./ rowSums(~zeroRows);
                probabilities(zeroRows, :) = 1 / nChoices;
            end
            estimator = struct( ...
                "kind", "categorical", ...
                "weights", ...
                    radia.optuna.internal.ParzenEstimator.mixtureWeights( ...
                    n, options.PriorWeight), ...
                "probabilities", probabilities, ...
                "n_choices", nChoices);
        end

        function values = sampleNumerical(estimator, stream, count)
            arguments
                estimator (1,1) struct
                stream (1,1) RandStream
                count (1,1) double {mustBeInteger, mustBePositive}
            end
            kernels = radia.optuna.internal.ParzenEstimator.sampleComponents( ...
                estimator.weights, stream, count);
            mu = estimator.mu(kernels);
            sigma = estimator.sigma(kernels);
            a = (estimator.internal_low - mu) ./ sigma;
            b = (estimator.internal_high - mu) ./ sigma;
            u = rand(stream, count, 1);
            cdfLow = radia.optuna.internal.ParzenEstimator.normalCdf(a);
            cdfHigh = radia.optuna.internal.ParzenEstimator.normalCdf(b);
            probabilities = cdfLow + u .* (cdfHigh - cdfLow);
            z = radia.optuna.internal.ParzenEstimator.normalInverse( ...
                probabilities);
            values = mu + sigma .* z;
            if estimator.log
                values = exp(values);
            end
            if isfinite(estimator.step)
                values = estimator.low + ...
                    round((values - estimator.low) / estimator.step) * ...
                    estimator.step;
            end
            values = min(max(values, estimator.low), estimator.high);
        end

        function values = sampleNumericalComponents(estimator, stream, components)
            % Sample specified mixture components; used by joint TPE.
            components = reshape(double(components), [], 1);
            if any(components < 1 | components > numel(estimator.weights) | ...
                    components ~= floor(components))
                error("radia:optuna:ParzenComponents", ...
                    "Invalid Parzen mixture component index.");
            end
            mu = estimator.mu(components);
            sigma = estimator.sigma(components);
            a = (estimator.internal_low - mu) ./ sigma;
            b = (estimator.internal_high - mu) ./ sigma;
            u = rand(stream, numel(components), 1);
            probabilities = radia.optuna.internal.ParzenEstimator.normalCdf(a) + ...
                u .* (radia.optuna.internal.ParzenEstimator.normalCdf(b) - ...
                radia.optuna.internal.ParzenEstimator.normalCdf(a));
            values = mu + sigma .* ...
                radia.optuna.internal.ParzenEstimator.normalInverse(probabilities);
            if estimator.log
                values = exp(values);
            end
            if isfinite(estimator.step)
                values = estimator.low + round((values - estimator.low) / ...
                    estimator.step) * estimator.step;
            end
            values = min(max(values, estimator.low), estimator.high);
        end

        function values = sampleCategorical(estimator, stream, count)
            arguments
                estimator (1,1) struct
                stream (1,1) RandStream
                count (1,1) double {mustBeInteger, mustBePositive}
            end
            kernels = radia.optuna.internal.ParzenEstimator.sampleComponents( ...
                estimator.weights, stream, count);
            values = zeros(count, 1);
            u = rand(stream, count, 1);
            for index = 1:count
                cumulative = cumsum(estimator.probabilities(kernels(index), :));
                cumulative(end) = 1;
                values(index) = 1 + sum(cumulative < u(index));
            end
        end

        function value = logPdfNumerical(estimator, samples)
            arguments
                estimator (1,1) struct
                samples double
            end
            samples = reshape(double(samples), [], 1);
            nSamples = numel(samples);
            nKernels = numel(estimator.weights);
            mu = reshape(estimator.mu, 1, nKernels);
            sigma = reshape(estimator.sigma, 1, nKernels);
            denominator = radia.optuna.internal.ParzenEstimator.logNormalMass( ...
                (estimator.internal_low - mu) ./ sigma, ...
                (estimator.internal_high - mu) ./ sigma);

            if isfinite(estimator.step)
                sampleLow = samples - estimator.step / 2;
                sampleHigh = samples + estimator.step / 2;
                if estimator.log
                    sampleLow = log(sampleLow);
                    sampleHigh = log(sampleHigh);
                end
                lower = (sampleLow - mu) ./ sigma;
                upper = (sampleHigh - mu) ./ sigma;
                componentLogPdf = ...
                    radia.optuna.internal.ParzenEstimator.logNormalMass( ...
                    lower, upper) - denominator;
            else
                internalSamples = samples;
                if estimator.log
                    internalSamples = log(samples);
                end
                z = (internalSamples - mu) ./ sigma;
                componentLogPdf = -0.5 * z.^2 - 0.5 * log(2*pi) - ...
                    log(sigma) - denominator;
            end
            if nSamples == 0
                value = zeros(0, 1);
                return
            end
            weighted = componentLogPdf + ...
                repmat(log(reshape(estimator.weights, 1, nKernels)), ...
                nSamples, 1);
            value = radia.optuna.internal.ParzenEstimator.logSumExp(weighted);
        end

        function value = logPdfCategorical(estimator, samples)
            arguments
                estimator (1,1) struct
                samples double
            end
            samples = reshape(double(samples), [], 1);
            if any(samples ~= floor(samples)) || ...
                    any(samples < 1 | samples > estimator.n_choices)
                error("radia:optuna:ParzenCategory", ...
                    "Categorical samples must index the available choices.");
            end
            nSamples = numel(samples);
            nKernels = numel(estimator.weights);
            probabilities = zeros(nSamples, nKernels);
            for index = 1:nSamples
                probabilities(index, :) = ...
                    estimator.probabilities(:, samples(index)).';
            end
            weighted = log(max(probabilities, realmin)) + ...
                repmat(log(reshape(estimator.weights, 1, nKernels)), ...
                nSamples, 1);
            value = radia.optuna.internal.ParzenEstimator.logSumExp(weighted);
        end

        function weights = defaultWeights(count)
            arguments
                count (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            if count == 0
                weights = zeros(0, 1);
            elseif count < 25
                weights = ones(count, 1);
            else
                rampCount = count - 25;
                if rampCount == 0
                    ramp = zeros(0, 1);
                else
                    ramp = linspace(1 / count, 1, rampCount).';
                end
                weights = [ramp; ones(25, 1)];
            end
        end
    end

    methods (Static, Access=private)
        function sigmas = sigmas(mus, low, high, magicClip, endpoints)
            count = numel(mus);
            if count == 0
                sigmas = zeros(0, 1);
                return
            end
            [sortedMus, order] = sort(mus);
            withEndpoints = [low; sortedMus; high];
            sortedSigmas = max( ...
                withEndpoints(2:end-1) - withEndpoints(1:end-2), ...
                withEndpoints(3:end) - withEndpoints(2:end-1));
            if ~endpoints && count >= 2
                sortedSigmas(1) = sortedMus(2) - sortedMus(1);
                sortedSigmas(end) = sortedMus(end) - sortedMus(end-1);
            end
            span = high - low;
            minimum = eps;
            if magicClip
                nKernels = count + 1;
                minimum = span / min(100, 1 + nKernels);
            end
            sortedSigmas = min(max(sortedSigmas, minimum), span);
            sigmas = zeros(count, 1);
            sigmas(order) = sortedSigmas;
        end

        function weights = mixtureWeights(count, priorWeight)
            if count == 0
                weights = 1;
                return
            end
            weights = [ ...
                radia.optuna.internal.ParzenEstimator.defaultWeights(count); ...
                priorWeight];
            total = sum(weights);
            if total <= 0
                error("radia:optuna:TPEWeights", ...
                    "Parzen mixture weights must contain positive mass.");
            end
            weights = weights / total;
        end

        function indices = sampleComponents(weights, stream, count)
            cumulative = cumsum(weights(:));
            cumulative(end) = 1;
            u = rand(stream, count, 1);
            indices = zeros(count, 1);
            for index = 1:count
                indices(index) = 1 + sum(cumulative < u(index));
            end
        end

        function values = normalCdf(z)
            values = 0.5 * erfc(-z / sqrt(2));
        end

        function values = normalInverse(probabilities)
            probabilities = min(max(probabilities, eps), 1 - eps);
            values = zeros(size(probabilities));
            lower = probabilities < 0.5;
            values(lower) = -sqrt(2) * erfcinv(2 * probabilities(lower));
            values(~lower) = sqrt(2) * ...
                erfcinv(2 * (1 - probabilities(~lower)));
        end

        function values = logNormalMass(a, b)
            mass = zeros(size(a));
            left = b <= 0;
            right = a >= 0;
            middle = ~(left | right);
            mass(left) = 0.5 * (erfc(-b(left) / sqrt(2)) - ...
                erfc(-a(left) / sqrt(2)));
            mass(right) = 0.5 * (erfc(a(right) / sqrt(2)) - ...
                erfc(b(right) / sqrt(2)));
            mass(middle) = 0.5 * (erf(b(middle) / sqrt(2)) - ...
                erf(a(middle) / sqrt(2)));
            values = log(max(mass, realmin));
        end

        function values = logSumExp(matrix)
            maximum = max(matrix, [], 2);
            finiteMaximum = maximum;
            finiteMaximum(~isfinite(finiteMaximum)) = 0;
            values = finiteMaximum + log(sum(exp(matrix - finiteMaximum), 2));
        end
    end
end
