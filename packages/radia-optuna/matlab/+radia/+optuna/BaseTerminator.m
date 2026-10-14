classdef (Abstract) BaseTerminator < handle
    %BASETERMINATOR Upstream Optuna automatic-stop interface.

    methods (Abstract)
        decision=should_terminate(obj,study)
    end
end
