classdef RegretBoundEvaluator < radia.optuna.BaseImprovementEvaluator
    %REGRETBOUNDEVALUATOR Gaussian-process regret bound from Optuna 4.9.

    properties (Access=private)
        PythonEvaluator
    end

    methods
        function obj=RegretBoundEvaluator(options)
            arguments
                options.top_trials_ratio (1,1) double {mustBePositive} = 0.5
                options.min_n_trials (1,1) double ...
                    {mustBeInteger,mustBePositive} = 20
                options.seed double = NaN
            end
            warning("radia:optuna:ExperimentalWarning", ...
                "RegretBoundEvaluator is experimental in Optuna 4.9.0.");
            module=py.importlib.import_module("optuna.terminator");
            evaluatorClass=py.builtins.getattr( ...
                module,"RegretBoundEvaluator");
            keyword={"top_trials_ratio",options.top_trials_ratio, ...
                "min_n_trials",int64(options.min_n_trials)};
            if ~isnan(options.seed)
                keyword(end+1:end+2)={"seed",int64(options.seed)};
            end
            obj.PythonEvaluator=evaluatorClass(pyargs(keyword{:}));
        end

        function value=evaluate(obj,trials,study_direction)
            [~,pythonTrials]=radia.optuna.internal.toUpstreamStudy( ...
                trials,study_direction);
            direction=obj.pythonDirection(study_direction);
            value=double(obj.PythonEvaluator.evaluate( ...
                pythonTrials,direction));
        end
    end

    methods (Static, Access=private)
        function direction=pythonDirection(value)
            module=py.importlib.import_module("optuna.study");
            directionClass=py.builtins.getattr(module,"StudyDirection");
            name=upper(string(radia.optuna.StudyDirection.from(value)));
            direction=py.builtins.getattr(directionClass,char(name));
        end
    end
end
