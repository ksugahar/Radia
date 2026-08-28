classdef LogUniformDistribution < radia.optuna.FloatDistribution
    %LOGUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.

    methods
        function obj=LogUniformDistribution(low,high)
            warning("radia:optuna:FutureWarning", ...
                "LogUniformDistribution is deprecated; use FloatDistribution(Log=true).");
            obj@radia.optuna.FloatDistribution(low,high,Log=true);
            obj.name="LogUniformDistribution";
        end
    end
end
