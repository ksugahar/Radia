classdef TrialState < uint8
    %TRIALSTATE Optuna 4.9 trial-state enumeration.

    enumeration
        RUNNING (0)
        COMPLETE (1)
        PRUNED (2)
        FAIL (3)
        WAITING (4)
    end

    methods
        function result=is_finished(obj)
            result=ismember(obj,[radia.optuna.TrialState.COMPLETE, ...
                radia.optuna.TrialState.PRUNED,radia.optuna.TrialState.FAIL]);
        end

        function text=string(obj)
            labels=["RUNNING","COMPLETE","PRUNED","FAIL","WAITING"];
            text=reshape(labels(double(obj)+1),size(obj));
        end

        function text=char(obj)
            if ~isscalar(obj)
                error("radia:optuna:TrialState", ...
                    "char requires a scalar TrialState.");
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
            if isa(value,"radia.optuna.TrialState")
                return
            end
            if isnumeric(value) && isscalar(value)
                members=[radia.optuna.TrialState.RUNNING, ...
                    radia.optuna.TrialState.COMPLETE, ...
                    radia.optuna.TrialState.PRUNED, ...
                    radia.optuna.TrialState.FAIL, ...
                    radia.optuna.TrialState.WAITING];
                match=find(double(members)==double(value),1);
                if isempty(match)
                    error("radia:optuna:TrialState", ...
                        "Unknown trial state.");
                end
                value=members(match);
                return
            end
            text=upper(string(value));
            valid=["RUNNING","COMPLETE","PRUNED","FAIL","WAITING"];
            if any(~ismember(text,valid))
                error("radia:optuna:TrialState","Unknown trial state.");
            end
            value=repmat(radia.optuna.TrialState.RUNNING,size(text));
            value(text=="COMPLETE")=radia.optuna.TrialState.COMPLETE;
            value(text=="PRUNED")=radia.optuna.TrialState.PRUNED;
            value(text=="FAIL")=radia.optuna.TrialState.FAIL;
            value(text=="WAITING")=radia.optuna.TrialState.WAITING;
        end

        function text=toStorage(value)
            text=string(radia.optuna.TrialState.from(value));
        end

        function value=from_bytes(bytes,byteorder,options)
            arguments
                bytes
                byteorder {mustBeTextScalar}
                options.Signed (1,1) logical = false
            end
            numeric=radia.optuna.internal.IntEnumSupport.fromBytes( ...
                bytes,byteorder,options.Signed);
            value=radia.optuna.TrialState.from(numeric);
        end
    end
end
