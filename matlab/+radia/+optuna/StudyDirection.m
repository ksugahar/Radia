classdef StudyDirection < uint8
    %STUDYDIRECTION Optuna 4.9 study-direction enumeration.

    enumeration
        NOT_SET (0)
        MINIMIZE (1)
        MAXIMIZE (2)
    end

    methods
        function text=string(obj)
            labels=["NOT_SET","MINIMIZE","MAXIMIZE"];
            text=reshape(labels(double(obj)+1),size(obj));
        end

        function text=char(obj)
            if ~isscalar(obj)
                error("radia:optuna:StudyDirection", ...
                    "char requires a scalar StudyDirection.");
            end
            text=char(string(obj));
        end
    end

    methods (Static)
        function value=from(value)
            if isa(value,"radia.optuna.StudyDirection")
                return
            end
            text=upper(string(value));
            if any(~ismember(text,["NOT_SET","MINIMIZE","MAXIMIZE"]))
                error("radia:optuna:Direction", ...
                    "Unknown study direction.");
            end
            value=repmat(radia.optuna.StudyDirection.NOT_SET,size(text));
            value(text=="MINIMIZE")=radia.optuna.StudyDirection.MINIMIZE;
            value(text=="MAXIMIZE")=radia.optuna.StudyDirection.MAXIMIZE;
        end

        function text=toStorage(value)
            value=radia.optuna.StudyDirection.from(value);
            text=lower(string(value));
            if any(text=="not_set")
                error("radia:optuna:Direction", ...
                    "NOT_SET is not a valid direction for create_study.");
            end
        end
    end
end
