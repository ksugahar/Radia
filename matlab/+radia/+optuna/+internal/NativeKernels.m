classdef NativeKernels
    %NATIVEKERNELS Enter the required, Optuna-only native gateway.
    methods (Static)
        function varargout = call(varargin)
            %CALL Enter the required Optuna-only MEX boundary.
            gateway = radia.optuna.internal.NativeKernels.gateway();
            if gateway == "optuna_mex"
                if ~isempty(varargin) && isstring(varargin{1}) && ...
                        isscalar(varargin{1})
                    varargin{1}=char(varargin{1});
                end
                [varargout{1:nargout}]=optuna_mex(varargin{:});
                return
            end
            error("radia:optuna:MissingNativeKernel", ...
                "optuna_mex is required. Run Build.ps1 -OptunaMexOnly " + ...
                "before using radia.optuna.");
        end

        function available = has(command)
            commands = radia.optuna.internal.NativeKernels.commands();
            available = any(commands == string(command));
        end

        function info = status()
            commands = radia.optuna.internal.NativeKernels.commands();
            expected = radia.optuna.internal.NativeKernels.expectedCommands();
            present = ismember(expected, commands);
            gateway=radia.optuna.internal.NativeKernels.gateway();
            if gateway == ""
                backend = "matlab-reference";
                mexPath = "";
            else
                backend = "native-mex";
                mexPath = string(which(gateway));
            end
            info = struct( ...
                schema="radia.optuna.native-kernels.v1", ...
                mex_available=~isempty(commands), ...
                backend=backend, ...
                gateway=gateway, ...
                mex_path=mexPath, ...
                expected_commands=expected, ...
                available_commands=expected(present), ...
                missing_commands=expected(~present));
        end
    end

    methods (Static, Access=private)
        function commands = commands()
            % A missing gateway reports an EMPTY command set rather than
            % raising, because every caller is written as
            %   if NativeKernels.has(cmd) ... else <MATLAB reference> end
            % and a raising predicate turns each of those guards into a hard
            % error while leaving the reference branch permanently
            % unreachable. The absence is announced once so the substitution
            % is never silent. A gateway that IS present but does not match
            % the checked command contract stays a hard error: a stale MEX
            % would compute silently different numbers.
            persistent cache initialized
            if isempty(initialized)
                gateway=radia.optuna.internal.NativeKernels.gateway();
                if gateway==""
                    radia.optuna.internal.NativeKernels.warnMissingGateway();
                    cache=strings(0,1);
                    initialized=true;
                    commands=cache;
                    return
                end
                cache=string(radia.optuna.internal.NativeKernels.call( ...
                    "api.commands"));
                expected=["api.info","api.commands", ...
                    radia.optuna.internal.NativeKernels.expectedCommands()];
                missing=expected(~ismember(expected,cache));
                unexpected=cache(~ismember(cache,expected));
                if numel(cache)~=numel(expected) || ...
                        ~isempty(missing) || ~isempty(unexpected)
                    error("radia:optuna:IncompatibleNativeKernel", ...
                        "optuna_mex does not match the required 20-command " + ...
                        "contract. Missing: %s. Unexpected: %s. Rebuild " + ...
                        "with Build.ps1 -OptunaMexOnly.", ...
                        strjoin(missing,", "),strjoin(unexpected,", "));
                end
                initialized=true;
            end
            commands=cache;
        end

        function warnMissingGateway()
            persistent warned
            if ~isempty(warned)
                return
            end
            warned = true;
            warning("radia:optuna:NativeKernelUnavailable", ...
                "optuna_mex is not on the MATLAB path, so radia.optuna " + ...
                "evaluates its MATLAB reference kernels instead of the " + ...
                "native gateway. Run Build.ps1 -OptunaMexOnly for the " + ...
                "supported configuration.");
        end

        function expected = expectedCommands()
            expected = [ ...
                "optuna.pareto.rank_crowding", ...
                "optuna.parzen.log_pdf_numerical", ...
                "optuna.parzen.log_pdf_categorical", ...
                "optuna.tpe.best_numerical", ...
                "optuna.tpe.best_joint", ...
                "optuna.tpe.best_joint_observations", ...
                "optuna.tpe.best_numerical_observations", ...
                "optuna.tpe.history.reset", ...
                "optuna.tpe.history.append_complete", ...
                "optuna.tpe.best_grouped_history", ...
                "optuna.random_state.create", ...
                "optuna.random_state.rand", ...
                "optuna.random_state.randn", ...
                "optuna.random_state.randi", ...
                "optuna.random_state.randperm", ...
                "optuna.random_state.snapshot", ...
                "optuna.random_state.restore", ...
                "optuna.random_state.destroy"];
        end

        function name = gateway()
            persistent cache initialized
            if isempty(initialized)
                if exist("optuna_mex","file") == 3
                    cache="optuna_mex";
                else
                    cache="";
                end
                initialized=true;
            end
            name=cache;
        end
    end
end
