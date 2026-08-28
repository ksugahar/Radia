classdef UniformDistribution < radia.optuna.FloatDistribution
    %UNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.

    methods
        function obj=UniformDistribution(low,high)
            warning("radia:optuna:FutureWarning", ...
                "UniformDistribution is deprecated; use FloatDistribution.");
            obj@radia.optuna.FloatDistribution(low,high);
            obj.name="UniformDistribution";
        end
    end
end
