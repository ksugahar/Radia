classdef HCurlTopologyOperator < handle
    %HCURLTOPOLOGYOPERATOR Matrix-free reduced HCurl HACApK operator.

    properties (SetAccess=private)
        ModeCount double = 0
        ChargeCount double = 0
        SubtetCount double = 0
        ParentCount double = 0
        Mu double = 0
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj=HCurlTopologyOperator(nativeHandle,info)
            obj.NativeHandle=uint64(nativeHandle);
            obj.ModeCount=info.mode_count;
            obj.ChargeCount=info.charge_count;
            obj.SubtetCount=info.subtet_count;
            obj.ParentCount=info.parent_count;
            obj.Mu=info.mu;
        end
    end

    methods (Static)
        function obj=create(gram,chargeMaps,cellVertices,chargeHosts,hostParents,options)
            arguments
                gram (1,1) radia.HACApKChargeGram
                chargeMaps double {mustBeFinite}
                cellVertices double {mustBeFinite}
                chargeHosts (:,1) {mustBeInteger,mustBeNonnegative}
                hostParents (:,1) {mustBeInteger,mustBeNonnegative}
                options.Mu (1,1) double {mustBePositive,mustBeFinite}=4*pi*1e-7
            end
            nativeHandle=radia.internal.callMex( ...
                'hcurl.topopt.operator.create',gram.nativeHandle(), ...
                double(chargeMaps),double(cellVertices),int32(chargeHosts), ...
                int32(hostParents),options.Mu);
            info=radia.internal.callMex('hcurl.topopt.operator.info',nativeHandle);
            obj=radia.topopt.HCurlTopologyOperator(nativeHandle,info);
        end
    end

    methods
        function value=info(obj)
            obj.assertAlive();
            value=radia.internal.callMex('hcurl.topopt.operator.info',obj.NativeHandle);
        end

        function value=matvec(obj,x)
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.matvec',obj.NativeHandle,double(x));
        end

        function value=toDense(obj)
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.to_dense',obj.NativeHandle);
        end

        function value=directionalContractions(obj,cellVertexVelocities,left,right)
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.directional_contractions', ...
                obj.NativeHandle,double(cellVertexVelocities),double(left),double(right));
        end

        function value=activationMatvec(obj,activation,x,options)
            arguments
                obj (1,1) radia.topopt.HCurlTopologyOperator
                activation (:,1) double {mustBeBetween(activation,0,1)}
                x double
                options.Power (1,1) double {mustBeGreaterThanOrEqual(options.Power,1)}=1
            end
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.activation_matvec',obj.NativeHandle, ...
                activation,double(x),options.Power);
        end

        function value=activationToDense(obj,activation,options)
            arguments
                obj (1,1) radia.topopt.HCurlTopologyOperator
                activation (:,1) double {mustBeBetween(activation,0,1)}
                options.Power (1,1) double {mustBeGreaterThanOrEqual(options.Power,1)}=1
            end
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.activation_to_dense',obj.NativeHandle, ...
                activation,options.Power);
        end

        function value=activationContractions(obj,activation,left,right,options)
            arguments
                obj (1,1) radia.topopt.HCurlTopologyOperator
                activation (:,1) double {mustBeBetween(activation,0,1)}
                left double
                right double
                options.Power (1,1) double {mustBeGreaterThanOrEqual(options.Power,1)}=1
            end
            obj.assertAlive();
            value=radia.internal.callMex( ...
                'hcurl.topopt.operator.activation_contractions',obj.NativeHandle, ...
                activation,double(left),double(right),options.Power);
        end

        function handle=nativeHandle(obj)
            obj.assertAlive();
            handle=obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle~=0
                try
                    radia.internal.callMex( ...
                        'hcurl.topopt.operator.destroy',obj.NativeHandle);
                catch
                end
                obj.NativeHandle=uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle==0
                error("radia:topopt:HCurlOperatorDeleted", ...
                    "The native HCurl topology operator has been deleted.");
            end
        end
    end
end
