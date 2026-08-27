classdef IntLogUniformDistribution < radia.optuna.IntDistribution
    %INTLOGUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.

    methods
        function obj=IntLogUniformDistribution(low,high,step)
            arguments
                low (1,1) double
                high (1,1) double
                step (1,1) double = 1
            end
            warning("radia:optuna:FutureWarning", ...
                "IntLogUniformDistribution is deprecated; use IntDistribution(Log=true).");
            obj@radia.optuna.IntDistribution(low,high,Log=true,Step=step);
            obj.name="IntLogUniformDistribution";
        end
    end
end
