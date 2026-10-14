classdef HACApKChargeGramDerivative < handle
    %HACAPKCHARGEGRAMDERIVATIVE Persistent directional-derivative H-matrix.

    properties (SetAccess=private)
        NDOF double = 0
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = HACApKChargeGramDerivative(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.NDOF = info.n_dof;
        end
    end

    methods (Static)
        function obj = fromNativeHandle(nativeHandle)
            info = radia.internal.callMex( ...
                'hacapk.charge_gram_derivative.info', nativeHandle);
            obj = radia.HACApKChargeGramDerivative(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'hacapk.charge_gram_derivative.info', obj.NativeHandle);
        end

        function value = stats(obj)
            value = obj.info();
        end

        function value = entry(obj, i, j)
            arguments
                obj (1,1) radia.HACApKChargeGramDerivative
                i (1,1) double {mustBeInteger,mustBePositive}
                j (1,1) double {mustBeInteger,mustBePositive}
            end
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'hacapk.charge_gram_derivative.entry', ...
                obj.NativeHandle, i, j);
        end

        function y = matvecSym(obj, x)
            arguments
                obj (1,1) radia.HACApKChargeGramDerivative
                x (:,1) double
            end
            obj.assertAlive();
            y = radia.internal.callMex( ...
                'hacapk.charge_gram_derivative.matvec_sym', ...
                obj.NativeHandle, x);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'hacapk.charge_gram_derivative.destroy', ...
                        obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HACApKChargeGramDerivative:Deleted", ...
                    "The native derivative operator has been deleted.");
            end
        end
    end
end
