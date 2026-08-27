classdef FloatDistribution < radia.optuna.BaseDistribution
    %FLOATDISTRIBUTION Public Optuna 4.9 float distribution.

    methods
        function obj=FloatDistribution(low,high,options)
            arguments
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = NaN
            end
            spec=radia.optuna.internal.DistributionCodec.float( ...
                low,high,options.Log,options.Step);
            if spec.high~=high
                warning("radia:optuna:DistributionAdjusted", ...
                    "FloatDistribution high was adjusted from %g to %g to align with Step=%g.", ...
                    high,spec.high,options.Step);
            end
            obj.name="FloatDistribution";
            obj.kind=spec.kind;
            obj.low=spec.low;
            obj.high=spec.high;
            obj.log=spec.log;
            obj.step=spec.step;
            obj.choices=spec.choices;
        end
    end
end
