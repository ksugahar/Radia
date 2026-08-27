classdef Terminator < radia.optuna.BaseTerminator
    %TERMINATOR Stop when estimated improvement no longer exceeds error.

    properties (SetAccess=private)
        ImprovementEvaluator
        ErrorEvaluator
        MinNTrials (1,1) double
    end

    methods
        function obj=Terminator(options)
            arguments
                options.ImprovementEvaluator = ...
                    radia.optuna.BestValueStagnationEvaluator()
                options.ErrorEvaluator = []
                options.MinNTrials (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 20
            end
            obj.ImprovementEvaluator=options.ImprovementEvaluator;
            obj.ErrorEvaluator=options.ErrorEvaluator;
            obj.MinNTrials=options.MinNTrials;
        end

        function decision=shouldTerminate(obj,study)
            complete=study.get_trials("COMPLETE");
            if numel(complete)<obj.MinNTrials
                decision=false;
                return
            end
            trials=study.get_trials();
            improvement=obj.ImprovementEvaluator.evaluate( ...
                trials,study.direction());
            if isempty(obj.ErrorEvaluator)
                errorValue=0;
            elseif isa(obj.ErrorEvaluator,"function_handle")
                errorValue=obj.ErrorEvaluator(study);
            else
                errorValue=obj.ErrorEvaluator.evaluate( ...
                    trials,study.direction());
            end
            decision=improvement<errorValue;
        end

        function decision=should_terminate(obj,study)
            decision=obj.shouldTerminate(study);
        end
    end
end
