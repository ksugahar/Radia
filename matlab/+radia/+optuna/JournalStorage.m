classdef JournalStorage < radia.optuna.RDBStorage
    %JOURNALSTORAGE Operation-log storage compatible with Optuna 4.9.
    %   The checked MATLAB file backend owns append/read interoperability;
    %   pinned Optuna owns replay ordering and conflict resolution.

    properties (SetAccess=immutable)
        LogStorage
    end

    methods
        function obj=JournalStorage(log_storage)
            arguments
                log_storage (1,1) radia.optuna.BaseJournalBackend
            end
            journalModule=py.importlib.import_module( ...
                "optuna.storages.journal");
            if isa(log_storage,"radia.optuna.JournalFileStorage")
                backendClass=py.builtins.getattr( ...
                    journalModule,"JournalFileBackend");
                pythonBackend=backendClass(char(log_storage.FilePath));
            elseif isa(log_storage,"radia.optuna.JournalRedisBackend")
                pythonBackend=log_storage.pythonBackend();
            else
                error("radia:optuna:JournalBackend", ...
                    "Unsupported journal backend '%s'.",class(log_storage));
            end
            storageClass=py.builtins.getattr(journalModule,"JournalStorage");
            pythonStorage=storageClass(pythonBackend);
            obj@radia.optuna.RDBStorage( ...
                "journal://"+string(class(log_storage)), ...
                InternalPythonStorage=pythonStorage);
            obj.LogStorage=log_storage;
        end

        function restore_replay_result(obj,snapshot)
            arguments
                obj
                snapshot {mustBeNumericOrLogical}
            end
            bytes=py.bytes(py.list(num2cell(reshape(uint8(snapshot),1,[]))));
            obj.PythonStorage.restore_replay_result(bytes);
        end
    end
end
