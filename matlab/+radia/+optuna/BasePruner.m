classdef (Abstract) BasePruner < handle
    %BASEPRUNER Common Optuna 4.9 pruning interface.

    methods
        function decision=prune(obj,study,trial)
            decision=obj.shouldPrune(study,trial); %#ok<MCNPN>
        end
    end
end
