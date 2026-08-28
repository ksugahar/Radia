classdef OptunaError < radia.optuna.ExceptionValue
    %OPTUNAERROR Base error corresponding to optuna.exceptions.OptunaError.

    methods
        function obj=OptunaError(varargin)
            obj@radia.optuna.ExceptionValue( ...
                "radia:optuna:OptunaError",varargin{:});
        end
    end
end
