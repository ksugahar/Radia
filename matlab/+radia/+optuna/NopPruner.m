classdef NopPruner
    %NOPPRUNER Never prune a trial.

    methods
        function decision=shouldPrune(~,~,~)
            decision=false;
        end
    end
end
