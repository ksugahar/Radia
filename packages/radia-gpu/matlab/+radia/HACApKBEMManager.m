classdef HACApKBEMManager < handle
    %HACAPKBEMMANAGER Stateful HACApK scalar BEM matrix.

    properties (SetAccess = private)
        NDOF
    end

    properties (Access = private)
        NativeHandle = uint64(0)
    end

    methods
        function obj = HACApKBEMManager(coordinates, entries)
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.bem.create', double(coordinates), double(entries));
            info = obj.info();
            obj.NDOF = info.n_dof;
        end

        function ok = build(obj, options)
            arguments
                obj
                options.AcaEps (1,1) double = -1
                options.LeafSize (1,1) double = -1
                options.Eta (1,1) double = -1
                options.MaxRank (1,1) double = -1
                options.PrintLevel (1,1) double = 0
            end
            obj.assertAlive();
            ok = radia.internal.callMex( ...
                'hacapk.bem.build', obj.NativeHandle, options.AcaEps, ...
                options.LeafSize, options.Eta, options.MaxRank, ...
                options.PrintLevel);
        end

        function y = matvec(obj, x)
            obj.assertAlive();
            y = radia.internal.callMex( ...
                'hacapk.bem.matvec', obj.NativeHandle, double(x));
        end

        function info = info(obj)
            obj.assertAlive();
            info = radia.internal.callMex('hacapk.bem.info', obj.NativeHandle);
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('hacapk.bem.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access = private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HACApKBEMManager:Deleted", ...
                    "The native HACApK BEM manager has been deleted.");
            end
        end
    end
end
