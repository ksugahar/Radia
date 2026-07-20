classdef RadiaField < radia.ngsolve.CoefficientFunction
    %RADIAFIELD Radia source field as a native NGSolve CoefficientFunction.

    methods
        function obj = RadiaField(radiaObject, fieldType, options)
            arguments
                radiaObject (1,1) double {mustBeInteger, mustBePositive}
                fieldType (1,1) string = "b"
                options.Origin double = []
                options.UAxis double = []
                options.VAxis double = []
                options.WAxis double = []
                options.Precision double = []
                options.Units (1,1) string = "m"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.radia_field.create', radiaObject, char(fieldType), ...
                options.Origin, options.UAxis, options.VAxis, options.WAxis, ...
                options.Precision, char(options.Units));
            obj@radia.ngsolve.CoefficientFunction(nativeHandle);
        end

        function value = fieldInfo(obj)
            value = radia.internal.callMex( ...
                'ngsolve.radia_field.info', obj.nativeHandle());
        end

        function prepareCache(obj, points)
            arguments
                obj (1,1) radia.RadiaField
                points (:,3) double
            end
            radia.internal.callMex('ngsolve.radia_field.prepare_cache', ...
                obj.nativeHandle(), points);
        end

        function clearCache(obj)
            radia.internal.callMex('ngsolve.radia_field.clear_cache', ...
                obj.nativeHandle());
        end

        function value = cacheStats(obj)
            value = radia.internal.callMex( ...
                'ngsolve.radia_field.cache_stats', obj.nativeHandle());
        end

        function coefficient = asVoxelCoefficient(obj, mesh, resolution)
            arguments
                obj (1,1) radia.RadiaField
                mesh (1,1) radia.ngsolve.Mesh
                resolution (1,1) double {mustBeInteger, mustBePositive} = 61
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.radia_field.as_voxel_coefficient', ...
                obj.nativeHandle(), mesh.nativeHandle(), resolution);
            coefficient = radia.ngsolve.CoefficientFunction.fromNativeHandle( ...
                nativeHandle);
        end
    end
end
