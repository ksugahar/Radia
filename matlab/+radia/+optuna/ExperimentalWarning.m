classdef ExperimentalWarning < radia.optuna.ExceptionValue
    %EXPERIMENTALWARNING Warning category for experimental Optuna features.

    methods
        function obj=ExperimentalWarning(varargin)
            obj@radia.optuna.ExceptionValue( ...
                "radia:optuna:ExperimentalWarning",varargin{:});
        end

        function warn(obj)
            warning(char(obj.identifier),"%s",obj.message);
        end
    end
end
