classdef PedAnovaImportanceEvaluator < radia.optuna.BaseImportanceEvaluator
    %PEDANOVAIMPORTANCEEVALUATOR Optuna 5.0 PED-ANOVA evaluator.

    properties (SetAccess=private)
        TargetQuantile (1,1) double = 0.1
        RegionQuantile (1,1) double = 1.0
        EvaluateOnLocal (1,1) logical = true
    end

    methods
        function obj=PedAnovaImportanceEvaluator(options)
            arguments
                options.target_quantile (1,1) double = 0.1
                options.region_quantile (1,1) double = 1.0
                options.evaluate_on_local (1,1) logical = true
            end
            obj.TargetQuantile=options.target_quantile;
            obj.RegionQuantile=options.region_quantile;
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
                EvaluateOnLocal=obj.EvaluateOnLocal);
        end
    end
end
