classdef JournalFileOpenLock < handle
    %JOURNALFILEOPENLOCK Cross-process exclusive journal lock.

    properties (SetAccess=immutable)
        FilePath (1,1) string
        GracePeriod double
    end

    properties (Access=private)
        LockPath (1,1) string
        Owned (1,1) logical = false
    end

    methods
        function obj=JournalFileOpenLock(filepath,grace_period)
            arguments
                filepath (1,1) string
                grace_period double = 30
            end
            if ~isempty(grace_period) && ...
                    (~isscalar(grace_period) || grace_period<=0 || ...
                    grace_period~=floor(grace_period))
                error("radia:optuna:JournalLockGrace", ...
                    "grace_period must be a positive integer or empty.");
            end
            if ~isempty(grace_period) && grace_period<3
                warning("radia:optuna:OptunaWarning", ...
                    "The value of grace_period might be too small.");
            end
            obj.FilePath=filepath;
            obj.LockPath=filepath+".lock";
            obj.GracePeriod=grace_period;
        end

        function acquired=acquire(obj)
            delay=0.001;
            observedTime=NaN;
            unchangedSince=tic;
            while true
                lockFile=java.io.File(char(obj.LockPath));
                if lockFile.createNewFile()
                    obj.Owned=true;
                    acquired=true;
                    return
                end
                if ~isempty(obj.GracePeriod) && isfile(obj.LockPath)
                    information=dir(obj.LockPath);
                    currentTime=information.datenum;
                    if ~isequaln(currentTime,observedTime)
                        observedTime=currentTime;
                        unchangedSince=tic;
                    elseif toc(unchangedSince)>obj.GracePeriod
                        delete(obj.LockPath);
                        delay=0.001;
                        continue
                    end
                end
                pause(delay);
                delay=min(delay*2,1);
            end
        end

        function release(obj)
            if ~obj.Owned || ~isfile(obj.LockPath)
                error("radia:optuna:JournalLockOwnership", ...
                    "Error: did not possess lock");
            end
            renamed=obj.LockPath+"."+string(java.util.UUID.randomUUID())+ ...
                ".rename";
            [moved,message]=movefile(obj.LockPath,renamed,"f");
            if ~moved
                error("radia:optuna:JournalLockOwnership","%s",message);
            end
            delete(renamed);
            obj.Owned=false;
        end
    end
end
