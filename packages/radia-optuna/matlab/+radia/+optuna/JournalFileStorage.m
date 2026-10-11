classdef JournalFileStorage < radia.optuna.BaseJournalBackend
    %JOURNALFILESTORAGE Newline-delimited JSON journal backend.

    properties (SetAccess=immutable)
        FilePath (1,1) string
        Lock
    end

    methods
        function obj=JournalFileStorage(file_path,lock_obj)
            arguments
                file_path (1,1) string
                lock_obj = []
            end
            warning("radia:optuna:FutureWarning", ...
                "JournalFileStorage is deprecated; use JournalFileBackend.");
            obj.FilePath=file_path;
            if isempty(lock_obj)
                lock_obj=radia.optuna.JournalFileSymlinkLock(file_path);
            elseif ~isa(lock_obj,"radia.optuna.JournalFileOpenLock")
                error("radia:optuna:JournalLock", ...
                    "lock_obj must implement the journal file lock contract.");
            end
            obj.Lock=lock_obj;
            if ~isfile(file_path)
                [fileId,message]=fopen(file_path,"a");
                if fileId<0
                    error("radia:optuna:JournalFile","%s",message);
                end
                fclose(fileId);
            end
        end

        function append_logs(obj,logs)
            entries=obj.logCells(logs);
            obj.Lock.acquire();
            cleanup=onCleanup(@()obj.releaseLock());
            [fileId,message]=fopen(obj.FilePath,"a","n","UTF-8");
            if fileId<0
                error("radia:optuna:JournalFile","%s",message);
            end
            fileCleanup=onCleanup(@()fclose(fileId));
            for index=1:numel(entries)
                fprintf(fileId,"%s\n",jsonencode(entries{index}));
            end
            clear fileCleanup
            obj.Lock.release();
            clear cleanup
        end

        function logs=read_logs(obj,log_number_from)
            arguments
                obj
                log_number_from (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            if ~isfile(obj.FilePath)
                logs=cell(0,1);
                return
            end
            bytes=readBytes(obj.FilePath);
            if ~isempty(bytes) && bytes(end)~=10
                error("radia:optuna:JournalFormat","Invalid log format.");
            end
            text=splitlines(string(native2unicode(bytes',"UTF-8")));
            if ~isempty(text) && text(end)==""
                text(end)=[];
            end
            first=log_number_from+1;
            if first>numel(text)
                logs=cell(0,1);
                return
            end
            selected=text(first:end);
            logs=cell(numel(selected),1);
            for index=1:numel(selected)
                try
                    logs{index}=jsondecode(char(selected(index)));
                catch cause
                    error("radia:optuna:JournalFormat", ...
                        "Invalid journal JSON: %s",cause.message);
                end
            end
        end
    end

    methods (Access=private)
        function releaseLock(obj)
            try
                obj.Lock.release();
            catch
                % Preserve an active append error; release is best effort here.
            end
        end

        function entries=logCells(~,logs)
            if iscell(logs)
                entries=reshape(logs,[],1);
            elseif isstruct(logs)
                entries=num2cell(reshape(logs,[],1));
            else
                error("radia:optuna:JournalLogs", ...
                    "logs must be a struct array or cell array of structs.");
            end
            if any(~cellfun(@(value)isstruct(value)&&isscalar(value),entries))
                error("radia:optuna:JournalLogs", ...
                    "Every journal log must be a scalar struct.");
            end
        end
    end
end

function bytes=readBytes(path)
fileId=fopen(path,"r");
if fileId<0
    error("radia:optuna:JournalFile","Cannot read journal file '%s'.",path);
end
cleanup=onCleanup(@()fclose(fileId));
bytes=fread(fileId,Inf,"*uint8");
clear cleanup
end
