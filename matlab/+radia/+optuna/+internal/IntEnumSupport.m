classdef IntEnumSupport
    %INTENUMSUPPORT Python IntEnum-compatible numeric helpers.

    methods (Static)
        function ratio=asIntegerRatio(value)
            ratio=[double(value),1];
        end

        function count=bitCount(value)
            count=arrayfun(@(item)sum(bitget(uint64(item),1:64)),double(value));
            count=reshape(count,size(value));
        end

        function count=bitLength(value)
            numeric=double(value);
            count=zeros(size(numeric));
            positive=numeric>0;
            count(positive)=floor(log2(numeric(positive)))+1;
        end

        function bytes=toBytes(value,length,byteorder)
            if ~isscalar(value) || length<0 || length~=fix(length)
                error("radia:optuna:IntEnumBytes", ...
                    "to_bytes requires a scalar value and nonnegative length.");
            end
            numeric=uint64(value);
            bytes=zeros(1,length,"uint8");
            for index=1:length
                bytes(length-index+1)=uint8(bitand(numeric,uint64(255)));
                numeric=bitshift(numeric,-8);
            end
            if numeric~=0
                error("radia:optuna:IntEnumBytes", ...
                    "Integer does not fit in the requested byte length.");
            end
            if lower(string(byteorder))=="little"
                bytes=fliplr(bytes);
            elseif lower(string(byteorder))~="big"
                error("radia:optuna:IntEnumBytes", ...
                    "byteorder must be 'big' or 'little'.");
            end
        end

        function numeric=fromBytes(bytes,byteorder,signed)
            bytes=reshape(uint8(bytes),1,[]);
            order=lower(string(byteorder));
            if order=="little"
                bytes=fliplr(bytes);
            elseif order~="big"
                error("radia:optuna:IntEnumBytes", ...
                    "byteorder must be 'big' or 'little'.");
            end
            numeric=uint64(0);
            for byte=bytes
                numeric=bitshift(numeric,8)+uint64(byte);
            end
            if signed && ~isempty(bytes) && bitget(bytes(1),8)
                numeric=double(numeric)-2^(8*numel(bytes));
            else
                numeric=double(numeric);
            end
        end
    end
end
