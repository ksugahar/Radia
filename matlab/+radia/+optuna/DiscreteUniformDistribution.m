classdef DiscreteUniformDistribution < radia.optuna.FloatDistribution
    %DISCRETEUNIFORMDISTRIBUTION Deprecated Optuna 4.9 distribution alias.

    properties (Dependent)
        q
    end

    methods
        function obj=DiscreteUniformDistribution(low,high,q)
            warning("radia:optuna:FutureWarning", ...
                "DiscreteUniformDistribution is deprecated; use FloatDistribution(Step=q).");
            obj@radia.optuna.FloatDistribution(low,high,Step=q);
            obj.name="DiscreteUniformDistribution";
        end

        function value=get.q(obj)
            value=obj.step;
        end
    end
end
