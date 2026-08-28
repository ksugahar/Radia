classdef JournalFileBackend < radia.optuna.JournalFileStorage
    %JOURNALFILEBACKEND Optuna 4.9 journal file backend name.

    methods
        function obj=JournalFileBackend(file_path,lock_obj)
            arguments
                file_path (1,1) string
                lock_obj = []
            end
            warningState=warning("off","radia:optuna:FutureWarning");
            cleanup=onCleanup(@()warning(warningState));
            obj@radia.optuna.JournalFileStorage(file_path,lock_obj);
            clear cleanup
        end
    end
end
