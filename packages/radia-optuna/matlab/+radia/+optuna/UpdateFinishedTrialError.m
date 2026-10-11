classdef UpdateFinishedTrialError < radia.optuna.OptunaError
    %UPDATEFINISHEDTRIALERROR A finished trial cannot be updated.

    methods
        function obj=UpdateFinishedTrialError(varargin)
            obj@radia.optuna.OptunaError(varargin{:});
            obj.identifier="radia:optuna:UpdateFinishedTrialError";
        end
    end
end
