classdef LogUniformDistribution < radia.optuna.FloatDistribution
    %LOGUNIFORMDISTRIBUTION Deprecated distribution alias retained by Optuna 5.

    methods
        function obj=LogUniformDistribution(low,high)
            warning("radia:optuna:FutureWarning", ...
                "LogUniformDistribution is deprecated; use FloatDistribution(Log=true).");
            obj@radia.optuna.FloatDistribution(low,high,Log=true);
            obj.name="LogUniformDistribution";
        end
    end
end
