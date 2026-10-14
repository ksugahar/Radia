classdef (Abstract) BaseImprovementEvaluator < handle
    %BASEIMPROVEMENTEVALUATOR Upstream Optuna improvement interface.

    methods (Abstract)
        value=evaluate(obj,trials,study_direction)
    end
end
