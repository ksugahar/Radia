classdef (Abstract) BaseCrossover
    %BASECROSSOVER MATLAB counterpart of Optuna's NSGA-II crossover base.

    properties (SetAccess=protected)
        NParents (1,1) double = 2
        Name (1,1) string = "BaseCrossover"
    end

    properties (Dependent, SetAccess=private)
        n_parents
    end

    methods
        function value=get.n_parents(obj)
            value=obj.NParents;
        end

        function config = configuration(obj)
            config = struct( ...
                "class",class(obj), ...
                "name",obj.Name, ...
                "n_parents",obj.NParents);
        end
    end

    methods (Abstract)
        child = crossover(obj,parents,stream,study,bounds)
    end
end
