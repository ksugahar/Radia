classdef NopPruner < radia.optuna.BasePruner
    %NOPPRUNER Never prune a trial.

    methods
        function decision=shouldPrune(~,~,~)
            decision=false;
        end
    end
end
