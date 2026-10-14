classdef UniformDistribution < radia.optuna.FloatDistribution
    %UNIFORMDISTRIBUTION Deprecated distribution alias retained by Optuna 5.

    methods
        function obj=UniformDistribution(low,high)
            warning("radia:optuna:FutureWarning", ...
                "UniformDistribution is deprecated; use FloatDistribution.");
            obj@radia.optuna.FloatDistribution(low,high);
            obj.name="UniformDistribution";
        end
    end
end
