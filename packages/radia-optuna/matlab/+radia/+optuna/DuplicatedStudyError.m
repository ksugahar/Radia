classdef DuplicatedStudyError < radia.optuna.OptunaError
    %DUPLICATEDSTUDYERROR A study with the requested name already exists.

    methods
        function obj=DuplicatedStudyError(varargin)
            obj@radia.optuna.OptunaError(varargin{:});
            obj.identifier="radia:optuna:DuplicatedStudyError";
        end
    end
end
