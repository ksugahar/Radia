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

        function result=as_integer_ratio(obj)
            result=radia.optuna.internal.IntEnumSupport.asIntegerRatio(obj);
        end

        function result=bit_count(obj)
            result=radia.optuna.internal.IntEnumSupport.bitCount(obj);
        end

        function result=bit_length(obj)
            result=radia.optuna.internal.IntEnumSupport.bitLength(obj);
        end

        function result=conjugate(obj)
            result=double(obj);
        end

        function result=denominator(~)
            result=1;
        end

        function result=imag(obj)
            result=zeros(size(obj));
        end

        function result=is_integer(obj)
            result=true(size(obj));
        end

        function result=name(obj)
            result=string(obj);
        end

        function result=numerator(obj)
            result=double(obj);
        end

        function result=real(obj)
            result=double(obj);
        end

        function bytes=to_bytes(obj,length,byteorder,options)
            arguments
                obj
                length (1,1) double
                byteorder {mustBeTextScalar}
                options.Signed (1,1) logical = false %#ok<INUSA>
            end
            bytes=radia.optuna.internal.IntEnumSupport.toBytes( ...
                obj,length,byteorder);
        end

        function result=value(obj)
            result=double(obj);
        end
    end

    methods (Static)
        function value=from(value)
            if isa(value,"radia.optuna.StudyDirection")
                return
            end
            if isnumeric(value) && isscalar(value)
                members=[radia.optuna.StudyDirection.NOT_SET, ...
                    radia.optuna.StudyDirection.MINIMIZE, ...
                    radia.optuna.StudyDirection.MAXIMIZE];
                match=find(double(members)==double(value),1);
                if isempty(match)
                    error("radia:optuna:Direction", ...
                        "Unknown study direction.");
                end
                value=members(match);
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

        function value=from_bytes(bytes,byteorder,options)
            arguments
                bytes
                byteorder {mustBeTextScalar}
                options.Signed (1,1) logical = false
            end
            numeric=radia.optuna.internal.IntEnumSupport.fromBytes( ...
                bytes,byteorder,options.Signed);
            value=radia.optuna.StudyDirection.from(numeric);
        end
    end
end
