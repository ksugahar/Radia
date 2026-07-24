classdef CMAEvolutionStrategy < handle
    %CMAEVOLUTIONSTRATEGY Standard full-covariance CMA-ES state machine.

    properties (SetAccess=private)
        Dimension (1,1) double
        PopulationSize (1,1) double
        Generation (1,1) double = 0
        Mean (1,:) double
        Sigma (1,1) double
        Covariance double
        PSigma (1,:) double
        PC (1,:) double
        Bounds double
        Stream
    end

    properties (Access=private)
        Mu (1,1) double
        MuEff (1,1) double
        Weights (:,1) double
        CC (1,1) double
        C1 (1,1) double
        CMu (1,1) double
        CSigma (1,1) double
        DSigma (1,1) double
        ChiN (1,1) double
        B double = zeros(0,0)
        D (:,1) double = zeros(0,1)
        MaxResampling (1,1) double = 100
    end

    methods
        function obj = CMAEvolutionStrategy(mean, sigma, options)
            arguments
                mean (1,:) double
                sigma (1,1) double
                options.Bounds double = zeros(0,2)
                options.PopulationSize (1,1) double = 0
                options.Seed (1,1) double = 0
                options.Covariance double = zeros(0,0)
                options.MaxResampling (1,1) double = 100
            end
            if isempty(mean) || any(~isfinite(mean))
                error("radia:optuna:CMAMean", ...
                    "CMA-ES mean must be a finite nonempty row vector.");
            end
            if ~(isfinite(sigma) && sigma > 0)
                error("radia:optuna:CMASigma", ...
                    "CMA-ES sigma must be positive and finite.");
            end
            dimension = numel(mean);
            populationSize = options.PopulationSize;
            if populationSize == 0
                populationSize = 4 + floor(3 * log(dimension));
            end
            if populationSize < 2 || populationSize ~= floor(populationSize)
                error("radia:optuna:CMAPopulation", ...
                    "CMA-ES population size must be an integer of at least two.");
            end
            if options.MaxResampling < 0 || ...
                    options.MaxResampling ~= floor(options.MaxResampling)
                error("radia:optuna:CMAResampling", ...
                    "MaxResampling must be a nonnegative integer.");
            end

            bounds = options.Bounds;
            if isempty(bounds)
                bounds = [-inf(dimension,1), inf(dimension,1)];
            end
            if ~isequal(size(bounds), [dimension,2]) || ...
                    any(isnan(bounds), "all") || ...
                    any(bounds(:,1) > bounds(:,2)) || ...
                    any(mean(:) < bounds(:,1)) || ...
                    any(mean(:) > bounds(:,2))
                error("radia:optuna:CMABounds", ...
                    "CMA-ES bounds must contain the initial mean.");
            end

            covariance = options.Covariance;
            if isempty(covariance)
                covariance = eye(dimension);
            end
            if ~isequal(size(covariance), [dimension,dimension]) || ...
                    any(~isfinite(covariance), "all")
                error("radia:optuna:CMACovariance", ...
                    "CMA-ES covariance has an invalid shape or value.");
            end

            obj.Dimension = dimension;
            obj.PopulationSize = populationSize;
            obj.Mean = double(mean);
            obj.Sigma = double(sigma);
            obj.Covariance = (double(covariance) + double(covariance).') / 2;
            obj.PSigma = zeros(1, dimension);
            obj.PC = zeros(1, dimension);
            obj.Bounds = double(bounds);
            obj.Stream = RandStream("mt19937ar", "Seed", double(options.Seed));
            obj.MaxResampling = options.MaxResampling;
            obj.initializeCoefficients();
            obj.repairCovariance();
        end

        function point = ask(obj)
            for attempt = 0:obj.MaxResampling
                point = obj.samplePoint();
                if all(point(:) >= obj.Bounds(:,1)) && ...
                        all(point(:) <= obj.Bounds(:,2))
                    return
                end
            end
            point = min(max(point, obj.Bounds(:,1).'), obj.Bounds(:,2).');
        end

        function tell(obj, points, fitness)
            points = double(points);
            fitness = reshape(double(fitness), [], 1);
            if ~isequal(size(points), [obj.PopulationSize,obj.Dimension]) || ...
                    numel(fitness) ~= obj.PopulationSize || ...
                    any(~isfinite(points), "all") || any(~isfinite(fitness))
                error("radia:optuna:CMATell", ...
                    "CMA-ES tell requires one finite full population.");
            end
            [~, order] = sort(fitness, "ascend");
            points = points(order,:);
            y = (points - obj.Mean) / obj.Sigma;
            yWeighted = obj.Weights(1:obj.Mu).' * y(1:obj.Mu,:);

            [basis, scales] = obj.eigendecomposition();
            inverseSqrt = basis * diag(1 ./ scales) * basis.';
            obj.Mean = obj.Mean + obj.Sigma * yWeighted;
            obj.PSigma = (1 - obj.CSigma) * obj.PSigma + ...
                sqrt(obj.CSigma * (2 - obj.CSigma) * obj.MuEff) * ...
                (inverseSqrt * yWeighted.').';

            normPSigma = norm(obj.PSigma);
            obj.Sigma = obj.Sigma * exp((obj.CSigma / obj.DSigma) * ...
                (normPSigma / obj.ChiN - 1));
            obj.Sigma = min(obj.Sigma, 1e32);
            obj.Generation = obj.Generation + 1;

            normalizedPath = normPSigma / sqrt(max(eps, ...
                1 - (1 - obj.CSigma)^(2 * (obj.Generation + 1))));
            threshold = (1.4 + 2 / (obj.Dimension + 1)) * obj.ChiN;
            hSigma = double(normalizedPath < threshold);
            obj.PC = (1 - obj.CC) * obj.PC + ...
                hSigma * sqrt(obj.CC * (2 - obj.CC) * obj.MuEff) * ...
                yWeighted;

            transformed = inverseSqrt * y.';
            transformedNorm2 = sum(transformed.^2, 1).';
            adjustedWeights = obj.Weights;
            negative = adjustedWeights < 0;
            adjustedWeights(negative) = adjustedWeights(negative) .* ...
                (obj.Dimension ./ (transformedNorm2(negative) + 1e-8));
            rankMu = zeros(obj.Dimension);
            for index = 1:obj.PopulationSize
                rankMu = rankMu + adjustedWeights(index) * ...
                    (y(index,:).' * y(index,:));
            end
            deltaHSigma = (1 - hSigma) * obj.CC * (2 - obj.CC);
            scale = 1 + obj.C1 * deltaHSigma - obj.C1 - ...
                obj.CMu * sum(obj.Weights);
            obj.Covariance = scale * obj.Covariance + ...
                obj.C1 * (obj.PC.' * obj.PC) + obj.CMu * rankMu;
            obj.repairCovariance();
            obj.B = zeros(0,0);
            obj.D = zeros(0,1);
        end

        function state = snapshot(obj)
            state = struct( ...
                "schema", "radia.optuna.cma-evolution-state.v1", ...
                "dimension", obj.Dimension, ...
                "population_size", obj.PopulationSize, ...
                "generation", obj.Generation, ...
                "mean", obj.Mean, ...
                "sigma", obj.Sigma, ...
                "covariance", obj.Covariance, ...
                "p_sigma", obj.PSigma, ...
                "p_c", obj.PC, ...
                "bounds", obj.Bounds, ...
                "max_resampling", obj.MaxResampling, ...
                "random_state", obj.Stream.State);
        end

        function result = shouldStop(obj)
            [~, scales] = obj.eigendecomposition();
            axisScale = obj.Sigma * scales;
            result = any(~isfinite(axisScale)) || ...
                max(axisScale) > 1e32 || ...
                max(scales) / max(min(scales), eps) > 1e14 || ...
                all(axisScale < 1e-12);
        end
    end

    methods (Static)
        function obj = fromSnapshot(state)
            required = ["schema","dimension","population_size", ...
                "generation","mean","sigma","covariance","p_sigma", ...
                "p_c","bounds","max_resampling","random_state"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state, required)) || ...
                    string(state.schema) ~= ...
                    "radia.optuna.cma-evolution-state.v1"
                error("radia:optuna:CMAState", ...
                    "CMA-ES state snapshot is invalid or unsupported.");
            end
            obj = radia.optuna.internal.CMAEvolutionStrategy( ...
                reshape(double(state.mean), 1, []), double(state.sigma), ...
                Bounds=double(state.bounds), ...
                PopulationSize=double(state.population_size), ...
                Covariance=double(state.covariance), ...
                MaxResampling=double(state.max_resampling));
            if obj.Dimension ~= double(state.dimension)
                error("radia:optuna:CMAState", ...
                    "CMA-ES state dimension is inconsistent.");
            end
            obj.Generation = double(state.generation);
            obj.PSigma = reshape(double(state.p_sigma), 1, []);
            obj.PC = reshape(double(state.p_c), 1, []);
            obj.Stream.State = state.random_state;
        end
    end

    methods (Access=private)
        function initializeCoefficients(obj)
            lambda = obj.PopulationSize;
            obj.Mu = floor(lambda / 2);
            raw = log((lambda + 1) / 2) - log((1:lambda).');
            positive = raw(1:obj.Mu);
            negative = raw(obj.Mu+1:end);
            obj.MuEff = sum(positive)^2 / sum(positive.^2);
            muEffMinus = sum(negative)^2 / sum(negative.^2);
            obj.C1 = 2 / ((obj.Dimension + 1.3)^2 + obj.MuEff);
            obj.CMu = min(1 - obj.C1 - 1e-8, ...
                2 * (obj.MuEff - 2 + 1 / obj.MuEff) / ...
                ((obj.Dimension + 2)^2 + obj.MuEff));
            minAlpha = min([1 + obj.C1 / obj.CMu, ...
                1 + 2 * muEffMinus / (obj.MuEff + 2), ...
                (1 - obj.C1 - obj.CMu) / ...
                (obj.Dimension * obj.CMu)]);
            raw(raw >= 0) = raw(raw >= 0) / sum(raw(raw >= 0));
            raw(raw < 0) = minAlpha * raw(raw < 0) / ...
                sum(abs(raw(raw < 0)));
            obj.Weights = raw;
            obj.CSigma = (obj.MuEff + 2) / ...
                (obj.Dimension + obj.MuEff + 5);
            obj.DSigma = 1 + 2 * max(0, ...
                sqrt((obj.MuEff - 1) / (obj.Dimension + 1)) - 1) + ...
                obj.CSigma;
            obj.CC = (4 + obj.MuEff / obj.Dimension) / ...
                (obj.Dimension + 4 + 2 * obj.MuEff / obj.Dimension);
            obj.ChiN = sqrt(obj.Dimension) * ...
                (1 - 1 / (4 * obj.Dimension) + ...
                1 / (21 * obj.Dimension^2));
        end

        function point = samplePoint(obj)
            [basis, scales] = obj.eigendecomposition();
            z = randn(obj.Stream, obj.Dimension, 1);
            point = obj.Mean + ...
                obj.Sigma * (basis * (scales .* z)).';
        end

        function [basis, scales] = eigendecomposition(obj)
            if isempty(obj.B) || isempty(obj.D)
                covariance = (obj.Covariance + obj.Covariance.') / 2;
                [basis, eigenvalues] = eig(covariance, "vector");
                eigenvalues = max(real(eigenvalues), eps);
                obj.B = real(basis);
                obj.D = sqrt(eigenvalues);
                obj.Covariance = obj.B * diag(obj.D.^2) * obj.B.';
            end
            basis = obj.B;
            scales = obj.D;
        end

        function repairCovariance(obj)
            covariance = (obj.Covariance + obj.Covariance.') / 2;
            [basis, eigenvalues] = eig(covariance, "vector");
            eigenvalues = max(real(eigenvalues), eps);
            obj.Covariance = real(basis) * diag(eigenvalues) * real(basis).';
        end
    end
end
