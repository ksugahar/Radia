classdef UniformCrossover < radia.optuna.nsgaii.BaseCrossover
    %UNIFORMCROSSOVER Select every numerical gene from one of two parents.

    properties (SetAccess=private)
        SwappingProbability (1,1) double = 0.5
    end

    methods
        function obj = UniformCrossover(options)
            arguments
                options.SwappingProbability (1,1) double = 0.5
            end
            mustBeProbability(options.SwappingProbability, ...
                "SwappingProbability");
            obj.NParents = 2;
            obj.Name = "UniformCrossover";
            obj.SwappingProbability = options.SwappingProbability;
        end

        function child = crossover(obj,parents,stream,study,bounds) %#ok<INUSD>
            count = size(parents,2);
            selectSecond = rand(stream,1,count) >= ...
                obj.SwappingProbability;
            child = parents(1,:);
            child(selectSecond) = parents(2,selectSecond);
        end

        function config = configuration(obj)
            config = configuration@radia.optuna.nsgaii.BaseCrossover(obj);
            config.swapping_probability = obj.SwappingProbability;
        end
    end
end

function mustBeProbability(value,name)
if ~isfinite(value) || value < 0 || value > 1
    error("radia:optuna:NSGAIIProbability", ...
        "%s must be in [0,1].",name);
end
end
