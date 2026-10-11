classdef (Abstract) BaseErrorEvaluator < handle
    %BASEERROREVALUATOR Upstream Optuna terminator error interface.

    methods (Abstract)
        value=evaluate(obj,trials,study_direction)
    end
end
