classdef NativeKernels
    %NATIVEKERNELS Discover optional radia_mex optimizer acceleration.
    methods (Static)
        function available = has(command)
            commands = radia.optuna.internal.NativeKernels.commands();
            available = any(commands == string(command));
        end

        function info = status()
            commands = radia.optuna.internal.NativeKernels.commands();
            expected = [ ...
                "optuna.pareto.rank_crowding", ...
                "optuna.parzen.log_pdf_numerical", ...
                "optuna.parzen.log_pdf_categorical"];
            present = ismember(expected, commands);
            backend = "matlab-fallback";
            if all(present)
                backend = "native-mex";
            end
            info = struct( ...
                schema="radia.optuna.native-kernels.v1", ...
                mex_available=~isempty(commands), ...
                backend=backend, ...
                expected_commands=expected, ...
                available_commands=expected(present), ...
                missing_commands=expected(~present));
        end
    end

    methods (Static, Access=private)
        function commands = commands()
            persistent cache
            mexPath = string(which("radia_mex"));
            if mexPath == ""
                commands = strings(0, 1);
                return
            end
            info = dir(mexPath);
            signature = mexPath + "|" + string(info.bytes) + "|" + ...
                compose("%.17g", info.datenum);
            if isempty(cache) || cache.Signature ~= signature
                cache = struct( ...
                    Signature=signature, ...
                    Commands=string(radia.internal.callMex("api.commands")));
            end
            commands = cache.Commands;
        end
    end
end
