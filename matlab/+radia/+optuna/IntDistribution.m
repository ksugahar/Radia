classdef IntDistribution < radia.optuna.BaseDistribution
    %INTDISTRIBUTION Public Optuna 4.9 integer distribution.

    methods
        function obj=IntDistribution(low,high,options)
            arguments
                low (1,1) double
                high (1,1) double
                options.Log (1,1) logical = false
                options.Step (1,1) double = 1
            end
            if low~=floor(low) || high~=floor(high) || ...
                    options.Step~=floor(options.Step) || options.Step<=0
                error("radia:optuna:IntegerDistribution", ...
                    "Integer bounds and Step must be integers and Step must be positive.");
            end
            if options.Log && options.Step~=1
                error("radia:optuna:IntegerDistribution", ...
                    "Log integer distributions require Step=1.");
            end
            spec=radia.optuna.internal.DistributionCodec.integer( ...
                low,high,options.Log,options.Step);
            if spec.high~=high
                warning("radia:optuna:DistributionAdjusted", ...
                    "IntDistribution high was adjusted from %g to %g to align with Step=%g.", ...
                    high,spec.high,options.Step);
            end
            obj.name="IntDistribution";
            obj.kind=spec.kind;
            obj.low=spec.low;
            obj.high=spec.high;
            obj.log=spec.log;
            obj.step=spec.step;
            obj.choices=spec.choices;
        end
    end
end
