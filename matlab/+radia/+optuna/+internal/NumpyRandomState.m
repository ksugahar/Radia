classdef NumpyRandomState < handle
    %NUMPYRANDOMSTATE NumPy RandomState-compatible MT19937 uniform stream.
    %   Optuna 4.9 samplers use NumPy's legacy RandomState. MATLAB's
    %   mt19937ar has the same engine family but a different seeding and
    %   floating-point extraction contract, so it cannot be used for a
    %   seeded differential oracle.

    properties (Access=private)
        MT (1,624) uint32 = zeros(1,624,"uint32")
        Index (1,1) double = 624
        HasGauss (1,1) logical = false
        Gauss (1,1) double = 0
        NativeHandle (1,1) uint64 = uint64(0)
        UseNative (1,1) logical = false
    end

    properties (Dependent)
        State
    end

    methods
        function obj=NumpyRandomState(seed)
            arguments
                seed (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            if seed>double(intmax("uint32"))
                error("radia:optuna:Seed", ...
                    "Seed must be in the uint32 range.");
            end
            nativeCommands=["optuna.random_state.create", ...
                "optuna.random_state.rand", "optuna.random_state.randn", ...
                "optuna.random_state.randi", "optuna.random_state.randperm", ...
                "optuna.random_state.snapshot", ...
                "optuna.random_state.restore", ...
                "optuna.random_state.destroy"];
            if all(arrayfun(@(command) ...
                    radia.optuna.internal.NativeKernels.has(command), ...
                    nativeCommands))
                obj.NativeHandle=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.create",double(seed));
                obj.UseNative=true;
                return
            end
            obj.MT(1)=uint32(seed);
            mask=uint64(4294967295);
            multiplier=uint64(1812433253);
            for index=2:624
                previous=obj.MT(index-1);
                mixed=bitxor(previous,bitshift(previous,-30));
                value=multiplier*uint64(mixed)+uint64(index-1);
                obj.MT(index)=uint32(bitand(value,mask));
            end
        end

        function value=get.State(obj)
            if obj.UseNative
                value=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.snapshot",obj.NativeHandle);
                return
            end
            value=struct("schema","numpy.randomstate.mt19937.v1", ...
                "mt",obj.MT,"index",obj.Index, ...
                "has_gauss",obj.HasGauss,"gauss",obj.Gauss);
        end

        function set.State(obj,value)
            if ~isstruct(value) || ~isscalar(value) || ...
                    ~all(isfield(value,["schema","mt","index"])) || ...
                    string(value.schema)~="numpy.randomstate.mt19937.v1" || ...
                    numel(value.mt)~=624 || value.index<0 || value.index>624
                error("radia:optuna:RandomState", ...
                    "Invalid NumPy RandomState snapshot.");
            end
            if isfield(value,"has_gauss")
                hasGauss=logical(value.has_gauss);
                gauss=double(value.gauss);
            else
                hasGauss=false;
                gauss=0;
                value.has_gauss=hasGauss;
                value.gauss=gauss;
            end
            if obj.UseNative
                radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.restore",obj.NativeHandle,value);
                return
            end
            obj.MT=reshape(uint32(value.mt),1,624);
            obj.Index=double(value.index);
            obj.HasGauss=hasGauss;
            obj.Gauss=gauss;
        end

        function delete(obj)
            if ~obj.UseNative || obj.NativeHandle==0
                return
            end
            try
                radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.destroy",obj.NativeHandle);
            catch
                % The MEX exit handler already owns cleanup after clear/unload.
            end
            obj.NativeHandle=uint64(0);
            obj.UseNative=false;
        end

        function values=rand(obj,varargin)
            if isempty(varargin)
                shape=[1,1];
            elseif isscalar(varargin) && ~isscalar(varargin{1})
                shape=double(varargin{1});
            else
                shape=cellfun(@double,varargin);
            end
            if isempty(shape), shape=[1,1]; end
            if obj.UseNative
                outputShape=obj.outputShape(shape);
                values=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.rand",obj.NativeHandle, ...
                    double(prod(outputShape)));
                values=reshape(values,outputShape);
                return
            end
            values=zeros(shape);
            for index=1:numel(values)
                left=bitshift(obj.nextUInt32(),-5);
                right=bitshift(obj.nextUInt32(),-6);
                values(index)=(double(left)*67108864+double(right))/9007199254740992;
            end
        end

        function values=randn(obj,varargin)
            if isempty(varargin)
                shape=[1,1];
            elseif isscalar(varargin) && ~isscalar(varargin{1})
                shape=double(varargin{1});
            else
                shape=cellfun(@double,varargin);
            end
            if isempty(shape), shape=[1,1]; end
            if obj.UseNative
                outputShape=obj.outputShape(shape);
                values=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.randn",obj.NativeHandle, ...
                    double(prod(outputShape)));
                values=reshape(values,outputShape);
                return
            end
            values=zeros(shape);
            for index=1:numel(values)
                if obj.HasGauss
                    values(index)=obj.Gauss;
                    obj.HasGauss=false;
                    obj.Gauss=0;
                    continue
                end
                radiusSquared=2;
                while radiusSquared>=1 || radiusSquared==0
                    left=2*rand(obj)-1;
                    right=2*rand(obj)-1;
                    radiusSquared=left*left+right*right;
                end
                factor=sqrt(-2*log(radiusSquared)/radiusSquared);
                obj.Gauss=factor*left;
                obj.HasGauss=true;
                values(index)=factor*right;
            end
        end

        function values=randi(obj,maximum,varargin)
            arguments
                obj
                maximum (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Repeating)
                varargin
            end
            if isempty(varargin)
                shape=[1,1];
            elseif isscalar(varargin) && ~isscalar(varargin{1})
                shape=double(varargin{1});
            else
                shape=cellfun(@double,varargin);
            end
            if obj.UseNative
                outputShape=obj.outputShape(shape);
                values=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.randi",obj.NativeHandle, ...
                    double(maximum),double(prod(outputShape)));
                values=reshape(values,outputShape);
                return
            end
            values=zeros(shape);
            for index=1:numel(values)
                values(index)=double(obj.interval(uint32(maximum-1)))+1;
            end
        end

        function values=randperm(obj,count)
            arguments
                obj
                count (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            if obj.UseNative
                values=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.random_state.randperm",obj.NativeHandle, ...
                    double(count));
                return
            end
            values=1:count;
            for right=count:-1:2
                left=double(obj.interval(uint32(right-1)))+1;
                temporary=values(right);
                values(right)=values(left);
                values(left)=temporary;
            end
        end

        function handle=nativeHandle(obj)
            if obj.UseNative
                handle=obj.NativeHandle;
            else
                handle=uint64(0);
            end
        end
    end

    methods (Access=private)
        function shape=outputShape(~,shape)
            if isempty(shape)
                shape=[1,1];
            elseif isscalar(shape)
                shape=[shape,shape];
            end
        end

        function value=interval(obj,maximum)
            mask=maximum;
            mask=bitor(mask,bitshift(mask,-1));
            mask=bitor(mask,bitshift(mask,-2));
            mask=bitor(mask,bitshift(mask,-4));
            mask=bitor(mask,bitshift(mask,-8));
            mask=bitor(mask,bitshift(mask,-16));
            while true
                value=bitand(obj.nextUInt32(),mask);
                if value<=maximum
                    return
                end
            end
        end

        function value=nextUInt32(obj)
            if obj.Index>=624
                obj.twist();
            end
            value=obj.MT(obj.Index+1);
            obj.Index=obj.Index+1;
            value=bitxor(value,bitshift(value,-11));
            value=bitxor(value,bitand(bitshift(value,7),uint32(hex2dec("9D2C5680"))));
            value=bitxor(value,bitand(bitshift(value,15),uint32(hex2dec("EFC60000"))));
            value=bitxor(value,bitshift(value,-18));
        end

        function twist(obj)
            upper=uint32(hex2dec("80000000"));
            lower=uint32(hex2dec("7FFFFFFF"));
            matrixA=uint32(hex2dec("9908B0DF"));
            % NumPy's legacy RandomState uses the reference in-place
            % two-segment twist. The second segment deliberately reads
            % values already written by the first segment.
            for zeroIndex=0:(624-397-1)
                joined=bitor(bitand(obj.MT(zeroIndex+1),upper), ...
                    bitand(obj.MT(zeroIndex+2),lower));
                value=bitxor(obj.MT(zeroIndex+397+1), ...
                    bitshift(joined,-1));
                if bitand(joined,uint32(1))~=0
                    value=bitxor(value,matrixA);
                end
                obj.MT(zeroIndex+1)=value;
            end
            for zeroIndex=(624-397):(624-2)
                joined=bitor(bitand(obj.MT(zeroIndex+1),upper), ...
                    bitand(obj.MT(zeroIndex+2),lower));
                value=bitxor(obj.MT(zeroIndex+397-624+1), ...
                    bitshift(joined,-1));
                if bitand(joined,uint32(1))~=0
                    value=bitxor(value,matrixA);
                end
                obj.MT(zeroIndex+1)=value;
            end
            joined=bitor(bitand(obj.MT(624),upper), ...
                bitand(obj.MT(1),lower));
            value=bitxor(obj.MT(397),bitshift(joined,-1));
            if bitand(joined,uint32(1))~=0
                value=bitxor(value,matrixA);
            end
            obj.MT(624)=value;
            obj.Index=0;
        end
    end
end
