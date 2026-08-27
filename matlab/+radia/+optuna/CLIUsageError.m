classdef CLIUsageError < radia.optuna.OptunaError
    %CLIUSAGEERROR Invalid Optuna command-line usage.

    methods
        function obj=CLIUsageError(varargin)
            obj@radia.optuna.OptunaError(varargin{:});
            obj.identifier="radia:optuna:CLIUsageError";
        end
    end
end
