classdef JournalFileSymlinkLock < radia.optuna.JournalFileOpenLock
    %JOURNALFILESYMLINKLOCK Portable MATLAB journal lock.
    %   MATLAB uses the same atomic exclusive-create contract on Windows,
    %   where creating NFS-style symbolic links is privilege dependent.

    methods
        function obj=JournalFileSymlinkLock(filepath,grace_period)
            arguments
                filepath (1,1) string
                grace_period double = 30
            end
            obj@radia.optuna.JournalFileOpenLock(filepath,grace_period);
        end
    end
end
