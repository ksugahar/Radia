classdef CategoricalDistribution < radia.optuna.BaseDistribution
    %CATEGORICALDISTRIBUTION Public Optuna 4.9 categorical distribution.

    methods
        function obj=CategoricalDistribution(choices)
            spec=radia.optuna.internal.DistributionCodec.categorical(choices);
            obj.name="CategoricalDistribution";
            obj.kind=spec.kind;
            obj.low=spec.low;
            obj.high=spec.high;
            obj.log=spec.log;
            obj.step=spec.step;
            obj.choices=spec.choices;
        end
    end
end
