classdef GridFunction < handle
    %GRIDFUNCTION Native NGSolve GridFunction with MATLAB vector access.
    %   Mesh, finite-element space, element maps, and coefficient evaluation
    %   remain owned by NGSolve. MATLAB accesses the DoF vector explicitly.

    properties (SetAccess=private)
        Space string = ""
        Order double = 0
        DofCount double = 0
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = GridFunction(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Space = string(info.space);
            obj.Order = info.order;
            obj.DofCount = info.dof_count;
        end
    end

    methods (Static)
        function obj = create(volPath, space, order, options)
            arguments
                volPath (1,1) string
                space (1,1) string
                order (1,1) double {mustBeInteger, mustBePositive}
                options.NoGrads (1,1) logical = true
                options.Name (1,1) string = "gfu"
                options.Complex (1,1) logical = false
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.grid_function.create', char(volPath), char(space), ...
                order, options.NoGrads, char(options.Name), options.Complex);
            info = radia.internal.callMex( ...
                'ngsolve.grid_function.info', nativeHandle);
            obj = radia.ngsolve.GridFunction(nativeHandle, info);
        end

        function obj = fromFESpace(space, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                options.Name (1,1) string = "gfu"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.grid_function.from_fespace', space.nativeHandle(), ...
                char(options.Name));
            info = radia.internal.callMex( ...
                'ngsolve.grid_function.info', nativeHandle);
            obj = radia.ngsolve.GridFunction(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'ngsolve.grid_function.info', obj.NativeHandle);
        end

        function values = vector(obj)
            obj.assertAlive();
            values = radia.internal.callMex( ...
                'ngsolve.grid_function.vector', obj.NativeHandle);
        end

        function result = vectorHandle(obj, component)
            arguments
                obj (1,1) radia.ngsolve.GridFunction
                component (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            obj.assertAlive();
            result = radia.ngsolve.Vector.fromGridFunction(obj, component);
        end

        function setVector(obj, values)
            arguments
                obj (1,1) radia.ngsolve.GridFunction
                values double
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'ngsolve.grid_function.set_vector', obj.NativeHandle, values);
        end

        function interpolate(obj, coefficient)
            arguments
                obj (1,1) radia.ngsolve.GridFunction
                coefficient (1,1) radia.ngsolve.CoefficientFunction
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'ngsolve.grid_function.interpolate', obj.NativeHandle, ...
                coefficient.nativeHandle());
        end

        function coefficient = asCoefficient(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.grid_function.as_coefficient', obj.NativeHandle);
            coefficient = radia.ngsolve.CoefficientFunction.fromNativeHandle(nativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'ngsolve.grid_function.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:ngsolve:GridFunctionDeleted", ...
                    "The native NGSolve GridFunction has been deleted.");
            end
        end
    end
end
