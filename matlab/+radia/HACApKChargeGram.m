classdef HACApKChargeGram < handle
    %HACAPKCHARGEGRAM Stateful monopole HDiv charge Gram H-matrix.

    properties (SetAccess = private)
        NDOF = 0
    end

    properties (Access = private)
        NativeHandle = uint64(0)
    end

    methods
        function obj = HACApKChargeGram(centroids, measures, selfEnergy)
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_monopole', double(centroids), ...
                double(measures), double(selfEnergy));
        end

        function ok = build(obj, options)
            arguments
                obj
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
            end
            obj.assertAlive();
            ok = radia.internal.callMex( ...
                'hacapk.charge_gram.build', obj.NativeHandle, ...
                options.AcaEps, options.LeafSize, options.Eta, ...
                options.MaxRank, options.PrintLevel);
            if ok
                obj.NDOF = radia.internal.callMex( ...
                    'hacapk.charge_gram.info', obj.NativeHandle).n_dof;
            end
        end

        function y = matvec(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec', x);
        end

        function y = matvecTranspose(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec_transpose', x);
        end

        function y = matvecSym(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec_sym', x);
        end

        function value = entry(obj, i, j)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'hacapk.charge_gram.entry', obj.NativeHandle, i, j);
        end

        function info = info(obj)
            obj.assertAlive();
            info = radia.internal.callMex( ...
                'hacapk.charge_gram.info', obj.NativeHandle);
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'hacapk.charge_gram.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access = private)
        function y = apply(obj, command, x)
            obj.assertAlive();
            y = radia.internal.callMex(command, obj.NativeHandle, double(x));
        end

        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HACApKChargeGram:Deleted", ...
                    "The native HACApK charge-Gram manager has been deleted.");
            end
        end
    end
end
