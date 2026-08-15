classdef SPXCrossover < radia.optuna.nsgaii.BaseCrossover
    %SPXCROSSOVER Three-parent simplex crossover.

    properties (SetAccess=private)
        Epsilon (1,1) double = NaN
    end

    methods
        function obj = SPXCrossover(options)
            arguments
                options.Epsilon (1,1) double = NaN
            end
            obj.NParents = 3;
            obj.Name = "SPXCrossover";
            obj.Epsilon = options.Epsilon;
        end

        function child = crossover(obj,parents,stream,study,bounds) %#ok<INUSD>
            centroid = mean(parents,1);
            count = obj.NParents-1;
            powers = 1./(1:count);
            coefficients = rand(stream,1,count).^powers;
            epsilon = obj.Epsilon;
            if isnan(epsilon)
                epsilon = sqrt(size(bounds,1)+2);
            end
            expanded = centroid+epsilon*(parents-centroid);
            accumulator = 0;
            for index = 2:obj.NParents
                accumulator = coefficients(index-1)*( ...
                    expanded(index-1,:)-expanded(index,:)+accumulator);
            end
            child = expanded(end,:)+accumulator;
        end

        function config = configuration(obj)
            config = configuration@radia.optuna.nsgaii.BaseCrossover(obj);
            config.epsilon = obj.Epsilon;
        end
    end
end
