classdef BLXAlphaCrossover < radia.optuna.nsgaii.BaseCrossover
    %BLXALPHACROSSOVER Blend crossover matching Optuna's implementation.

    properties (SetAccess=private)
        Alpha (1,1) double = 0.5
    end

    methods
        function obj = BLXAlphaCrossover(options)
            arguments
                options.Alpha (1,1) double = 0.5
            end
            obj.NParents = 2;
            obj.Name = "BLXAlphaCrossover";
            obj.Alpha = options.Alpha;
        end

        function child = crossover(obj,parents,stream,study,bounds) %#ok<INUSD>
            parentsMin = min(parents,[],1);
            parentsMax = max(parents,[],1);
            extension = obj.Alpha*(parentsMax-parentsMin);
            low = parentsMin-extension;
            high = parentsMax+extension;
            child = low+(high-low).*rand(stream,1,size(parents,2));
        end

        function config = configuration(obj)
            config = configuration@radia.optuna.nsgaii.BaseCrossover(obj);
            config.alpha = obj.Alpha;
        end
    end
end
