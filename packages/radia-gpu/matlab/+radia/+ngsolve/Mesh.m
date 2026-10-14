classdef Mesh < handle
    %MESH Persistent native NGSolve MeshAccess handle.
    %   The mesh remains owned by NGSolve and is shared by dependent spaces.

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = Mesh(nativeHandle)
            obj.NativeHandle = uint64(nativeHandle);
        end
    end

    methods (Static)
        function obj = create(volPath)
            arguments
                volPath (1,1) string
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.mesh.create', char(volPath));
            obj = radia.ngsolve.Mesh(nativeHandle);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.mesh.info', ...
                obj.NativeHandle);
        end

        function setDeformation(obj, deformation)
            arguments
                obj (1,1) radia.ngsolve.Mesh
                deformation (1,1) radia.ngsolve.GridFunction
            end
            obj.assertAlive();
            radia.internal.callMex('ngsolve.mesh.set_deformation', ...
                obj.NativeHandle, deformation.nativeHandle());
        end

        function unsetDeformation(obj)
            obj.assertAlive();
            radia.internal.callMex('ngsolve.mesh.unset_deformation', ...
                obj.NativeHandle);
        end

        function quality = trafoQuality(obj, options)
            arguments
                obj (1,1) radia.ngsolve.Mesh
                options.IntegrationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 2
                options.ReferenceDeterminants (:,1) double = []
            end
            obj.assertAlive();
            if isempty(options.ReferenceDeterminants)
                quality = radia.internal.callMex( ...
                    'ngsolve.mesh.trafo_quality', obj.NativeHandle, ...
                    options.IntegrationOrder);
            else
                quality = radia.internal.callMex( ...
                    'ngsolve.mesh.trafo_quality', obj.NativeHandle, ...
                    options.IntegrationOrder, options.ReferenceDeterminants);
            end
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.mesh.destroy', ...
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
                error("radia:ngsolve:MeshDeleted", ...
                    "The native NGSolve Mesh has been deleted.");
            end
        end
    end
end
