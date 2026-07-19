classdef HACApKPEECManager < handle
    %HACAPKPEECMANAGER Stateful HACApK PEEC filament inductance matrix.

    properties (SetAccess = private)
        NDOF = 0
    end

    properties (Access = private)
        NativeHandle = uint64(0)
    end

    methods
        function obj = HACApKPEECManager(centers, directions, lengths, ...
                widths, heights, sigmas)
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.peec.create', double(centers), double(directions), ...
                double(lengths), double(widths), double(heights), double(sigmas));
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
                'hacapk.peec.build', obj.NativeHandle, options.AcaEps, ...
                options.LeafSize, options.Eta, options.MaxRank, ...
                options.PrintLevel);
            if ok
                obj.NDOF = radia.internal.callMex( ...
                    'hacapk.peec.info', obj.NativeHandle).n_dof;
            end
        end

        function y = matvec(obj, x)
            obj.assertAlive();
            y = radia.internal.callMex( ...
                'hacapk.peec.matvec', obj.NativeHandle, double(x));
        end

        function info = info(obj)
            obj.assertAlive();
            info = radia.internal.callMex('hacapk.peec.info', obj.NativeHandle);
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('hacapk.peec.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access = private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HACApKPEECManager:Deleted", ...
                    "The native HACApK PEEC manager has been deleted.");
            end
        end
    end
end
