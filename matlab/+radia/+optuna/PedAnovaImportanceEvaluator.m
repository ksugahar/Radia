classdef PedAnovaImportanceEvaluator < radia.optuna.BaseImportanceEvaluator
    %PEDANOVAIMPORTANCEEVALUATOR Optuna 4.9 PED-ANOVA evaluator.

    properties (SetAccess=private)
        TargetQuantile (1,1) double = 0.1
        RegionQuantile (1,1) double = 1.0
        BaselineQuantile (1,1) double = NaN
        EvaluateOnLocal (1,1) logical = true
    end

    methods
        function obj=PedAnovaImportanceEvaluator(options)
            arguments
                options.target_quantile (1,1) double = 0.1
                options.region_quantile (1,1) double = 1.0
                options.baseline_quantile (1,1) double = NaN
                options.evaluate_on_local (1,1) logical = true
            end
            obj.TargetQuantile=options.target_quantile;
            obj.RegionQuantile=options.region_quantile;
            obj.BaselineQuantile=options.baseline_quantile;
            obj.EvaluateOnLocal=options.evaluate_on_local;
        end

        function result=evaluate(obj,study,params,options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                params string = strings(1,0)
                options.target = []
            end
            result=obj.evaluateWith(study,"ped_anova",params,options.target, ...
                TargetQuantile=obj.TargetQuantile, ...
                RegionQuantile=obj.RegionQuantile, ...
                BaselineQuantile=obj.BaselineQuantile, ...
                EvaluateOnLocal=obj.EvaluateOnLocal);
        end
    end
end
