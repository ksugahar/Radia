classdef FESpace < handle
    %FESPACE Persistent native NGSolve finite-element space handle.

    properties (SetAccess=private)
        Space string = ""
        Order double = 0
        DofCount double = 0
        FreeDofCount double = 0
        Dirichlet string = ""
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = FESpace(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Space = string(info.space);
            obj.Order = info.order;
            obj.DofCount = info.dof_count;
            obj.FreeDofCount = info.free_dof_count;
            obj.Dirichlet = string(info.dirichlet);
        end
    end

    methods (Static)
        function obj = create(mesh, space, order, options)
            arguments
                mesh (1,1) radia.ngsolve.Mesh
                space (1,1) string
                order (1,1) double {mustBeInteger, mustBePositive}
                options.NoGrads (1,1) logical = true
                options.Complex (1,1) logical = false
                options.Dirichlet (1,1) string = ""
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.fespace.create', mesh.nativeHandle(), ...
                char(space), order, options.NoGrads, options.Complex, ...
                char(options.Dirichlet));
            info = radia.internal.callMex('ngsolve.fespace.info', ...
                nativeHandle);
            obj = radia.ngsolve.FESpace(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.fespace.info', ...
                obj.NativeHandle);
        end

        function value = freeDofs(obj)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'ngsolve.fespace.free_dofs', obj.NativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.fespace.destroy', ...
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
                error("radia:ngsolve:FESpaceDeleted", ...
                    "The native NGSolve FESpace has been deleted.");
            end
        end
    end
end
