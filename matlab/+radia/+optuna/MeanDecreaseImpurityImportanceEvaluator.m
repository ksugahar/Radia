classdef MeanDecreaseImpurityImportanceEvaluator < ...
        radia.optuna.BaseImportanceEvaluator
    %MEANDECREASEIMPURITYIMPORTANCEEVALUATOR Optuna 4.9 MDI evaluator.

    properties (SetAccess=private)
        NTrees (1,1) double = 64
        MaxDepth (1,1) double = 64
        Seed (1,1) double = 0
    end

    methods
        function obj=MeanDecreaseImpurityImportanceEvaluator(options)
            arguments
                options.n_trees (1,1) double = 64
                options.max_depth (1,1) double = 64
                options.seed (1,1) double = 0
            end
            obj.NTrees=options.n_trees;
            obj.MaxDepth=options.max_depth;
            obj.Seed=options.seed;
        end

        function result=evaluate(obj,study,params,options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                params string = strings(1,0)
                options.target = []
            end
            result=obj.evaluateWith(study,"mdi",params,options.target, ...
                Seed=obj.Seed,NTrees=obj.NTrees,MaxDepth=obj.MaxDepth);
        end
    end
end
