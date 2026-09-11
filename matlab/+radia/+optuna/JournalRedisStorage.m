classdef JournalRedisStorage < radia.optuna.JournalRedisBackend
    %JOURNALREDISSTORAGE Deprecated Redis backend name retained by Optuna 5.

    methods
        function obj=JournalRedisStorage(url,options)
            arguments
                url (1,1) string
                options.use_cluster (1,1) logical = false
                options.prefix (1,1) string = ""
                options.Client = []
            end
            warning("radia:optuna:FutureWarning", ...
                "JournalRedisStorage is deprecated; use JournalRedisBackend.");
            obj@radia.optuna.JournalRedisBackend(url, ...
                use_cluster=options.use_cluster,prefix=options.prefix, ...
                Client=options.Client);
        end
    end
end
