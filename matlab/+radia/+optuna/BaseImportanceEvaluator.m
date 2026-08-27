classdef (Abstract) BaseImportanceEvaluator
    %BASEIMPORTANCEEVALUATOR Public Optuna 4.9 importance interface.

    methods (Abstract)
        result=evaluate(obj,study,params,options)
    end

    methods (Access=protected)
        function result=evaluateWith(obj,study,evaluator,params,target,options) %#ok<INUSD>
            arguments
                obj
                study (1,1) radia.optuna.Study
                evaluator (1,1) string
                params string
                target
                options.Seed (1,1) double = 0
                options.NTrees (1,1) double = 64
                options.MaxDepth (1,1) double = 64
                options.TargetQuantile (1,1) double = 0.1
                options.RegionQuantile (1,1) double = 1.0
                options.BaselineQuantile (1,1) double = NaN
                options.EvaluateOnLocal (1,1) logical = true
            end
            result=radia.optuna.get_param_importances(study, ...
                evaluator=evaluator,params=params,normalize=false, ...
                seed=options.Seed,n_trees=options.NTrees, ...
                max_depth=options.MaxDepth, ...
                target_quantile=options.TargetQuantile, ...
                region_quantile=options.RegionQuantile, ...
                baseline_quantile=options.BaselineQuantile, ...
                evaluate_on_local=options.EvaluateOnLocal,target=target);
        end
    end
end
