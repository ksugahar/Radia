classdef IntUniformDistribution < radia.optuna.IntDistribution
    %INTUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.

    methods
        function obj=IntUniformDistribution(low,high,step)
            arguments
                low (1,1) double
                high (1,1) double
                step (1,1) double = 1
            end
            warning("radia:optuna:FutureWarning", ...
                "IntUniformDistribution is deprecated; use IntDistribution.");
            obj@radia.optuna.IntDistribution(low,high,Step=step);
            obj.name="IntUniformDistribution";
        end
    end
end
