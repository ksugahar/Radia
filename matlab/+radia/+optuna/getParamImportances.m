function result=getParamImportances(study,options)
%GETPARAMIMPORTANCES CamelCase alias for GET_PARAM_IMPORTANCES.
arguments
    study (1,1) radia.optuna.Study
    options.Evaluator (1,1) string = "fanova"
    options.Params string = strings(1,0)
    options.Normalize (1,1) logical = true
    options.ObjectiveIndex (1,1) double {mustBeInteger,mustBePositive} = 1
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.NTrees (1,1) double {mustBeInteger,mustBePositive} = 64
    options.MaxDepth (1,1) double {mustBeInteger,mustBePositive} = 64
    options.TargetQuantile (1,1) double = 0.1
    options.RegionQuantile (1,1) double = 1.0
end
result=radia.optuna.get_param_importances(study, ...
    evaluator=options.Evaluator,params=options.Params, ...
    normalize=options.Normalize,objective_index=options.ObjectiveIndex, ...
    seed=options.Seed,n_trees=options.NTrees,max_depth=options.MaxDepth, ...
    target_quantile=options.TargetQuantile, ...
    region_quantile=options.RegionQuantile);
end
