classdef StaticErrorEvaluator < radia.optuna.BaseErrorEvaluator
    %STATICERROREVALUATOR Return a fixed statistical-error estimate.

    properties (SetAccess=private)
        Constant (1,1) double
    end

    methods
        function obj=StaticErrorEvaluator(constant)
            arguments
                constant (1,1) double
            end
            obj.Constant=constant;
        end

        function value=evaluate(obj,trials,study_direction) %#ok<INUSD>
            value=obj.Constant;
        end
    end
end
