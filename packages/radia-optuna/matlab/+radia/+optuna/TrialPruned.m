classdef TrialPruned < radia.optuna.OptunaError
    %TRIALPRUNED Exception used to prune an objective trial.
    %   throw(radia.optuna.TrialPruned()) records a PRUNED trial, preserves
    %   the last reported intermediate value, and continues optimization.

    methods
        function obj=TrialPruned(varargin)
            obj@radia.optuna.OptunaError(varargin{:});
            obj.identifier="radia:optuna:TrialPruned";
        end
    end
end
