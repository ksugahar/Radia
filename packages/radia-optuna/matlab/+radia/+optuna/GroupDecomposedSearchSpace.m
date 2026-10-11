classdef GroupDecomposedSearchSpace < handle
    %GROUPDECOMPOSEDSEARCHSPACE MATLAB name for Optuna's underscored class.

    properties (SetAccess=private)
        IncludePruned (1,1) logical = false
    end

    properties (Access=private)
        Core (1,1) radia.optuna.internal.GroupDecomposedSearchSpace
    end

    methods
        function obj=GroupDecomposedSearchSpace(options)
            arguments
                options.include_pruned (1,1) logical = false
            end
            obj.IncludePruned=options.include_pruned;
            obj.Core=radia.optuna.internal.GroupDecomposedSearchSpace();
        end

        function result=calculate(obj,study,options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                options.use_cache (1,1) logical = false %#ok<INUSA>
            end
            groups=obj.Core.calculate(study, ...
                IncludePruned=obj.IncludePruned,ExcludeSingle=false);
            result=radia.optuna.SearchSpaceGroup(groups);
        end
    end
end
