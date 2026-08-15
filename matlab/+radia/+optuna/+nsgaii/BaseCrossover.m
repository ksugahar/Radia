classdef BaseCrossover
    %BASECROSSOVER MATLAB counterpart of Optuna's NSGA-II crossover base.

    properties (SetAccess=protected)
        NParents (1,1) double = 2
        Name (1,1) string = "BaseCrossover"
    end

    methods
        function child = crossover(~,parents,stream,study,bounds) %#ok<INUSD>
            child = zeros(1,0); %#ok<NASGU>
            error("radia:optuna:NSGAIICrossover", ...
                "Concrete crossover classes must implement crossover().");
        end

        function config = configuration(obj)
            config = struct( ...
                "class",class(obj), ...
                "name",obj.Name, ...
                "n_parents",obj.NParents);
        end
    end
end
