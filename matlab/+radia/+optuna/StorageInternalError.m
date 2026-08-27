classdef StorageInternalError < radia.optuna.OptunaError
    %STORAGEINTERNALERROR Optuna storage backend failure.

    methods
        function obj=StorageInternalError(varargin)
            obj@radia.optuna.OptunaError(varargin{:});
            obj.identifier="radia:optuna:StorageInternalError";
        end
    end
end
