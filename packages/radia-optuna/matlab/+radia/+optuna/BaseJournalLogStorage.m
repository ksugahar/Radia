classdef (Abstract) BaseJournalLogStorage < handle
    %BASEJOURNALLOGSTORAGE Deprecated Optuna journal backend interface.

    methods (Abstract)
        append_logs(obj,logs)
        logs=read_logs(obj,log_number_from)
    end
end
