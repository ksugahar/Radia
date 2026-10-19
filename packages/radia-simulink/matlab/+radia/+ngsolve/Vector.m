classdef Vector < handle
    %VECTOR Native NGSolve BaseVector for low-copy iterative workflows.
    %   A vector made from GridFunction.vectorHandle() is a live view of one
    %   GridFunction component. copy() creates an independent work vector;
    %   arithmetic then stays in the native NGSolve representation.

    properties (SetAccess=private)
        Size double = 0
        IsComplex logical = false
        IsView logical = false
        Component double = 0
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = Vector(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Size = info.size;
            obj.IsComplex = info.is_complex;
            obj.IsView = info.is_view;
            obj.Component = info.component;
        end
    end

    methods (Static)
        function obj = fromGridFunction(gridFunction, component)
            arguments
                gridFunction (1,1) radia.ngsolve.GridFunction
                component (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.grid_function.vector_handle', ...
                gridFunction.nativeHandle(), component);
            info = radia.internal.callMex('ngsolve.vector.info', nativeHandle);
            obj = radia.ngsolve.Vector(nativeHandle, info);
        end

        function obj = fromNativeHandle(nativeHandle)
            info = radia.internal.callMex('ngsolve.vector.info', nativeHandle);
            obj = radia.ngsolve.Vector(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.vector.info', ...
                obj.NativeHandle);
        end

        function result = copy(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex('ngsolve.vector.copy', ...
                obj.NativeHandle);
            result = radia.ngsolve.Vector.fromNativeHandle(nativeHandle);
        end

        function setZero(obj)
            obj.assertAlive();
            radia.internal.callMex('ngsolve.vector.set_zero', ...
                obj.NativeHandle);
        end

        function scale(obj, scalar)
            arguments
                obj (1,1) radia.ngsolve.Vector
                scalar (1,1) double
            end
            obj.assertAlive();
            radia.internal.callMex('ngsolve.vector.scale', scalar, ...
                obj.NativeHandle);
        end

        function axpy(obj, alpha, x)
            arguments
                obj (1,1) radia.ngsolve.Vector
                alpha (1,1) double
                x (1,1) radia.ngsolve.Vector
            end
            obj.assertAlive();
            x.assertAlive();
            radia.internal.callMex('ngsolve.vector.axpy', alpha, ...
                obj.NativeHandle, x.NativeHandle);
        end

        function value = dot(obj, other, options)
            arguments
                obj (1,1) radia.ngsolve.Vector
                other (1,1) radia.ngsolve.Vector
                options.Conjugate (1,1) logical = false
            end
            obj.assertAlive();
            other.assertAlive();
            value = radia.internal.callMex('ngsolve.vector.dot', ...
                obj.NativeHandle, other.NativeHandle, options.Conjugate);
        end

        function value = norm(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.vector.norm', ...
                obj.NativeHandle);
        end

        function values = values(obj)
            obj.assertAlive();
            values = radia.internal.callMex('ngsolve.vector.values', ...
                obj.NativeHandle);
        end

        function setValues(obj, values)
            arguments
                obj (1,1) radia.ngsolve.Vector
                values double
            end
            obj.assertAlive();
            radia.internal.callMex('ngsolve.vector.set_values', ...
                obj.NativeHandle, values);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.vector.destroy', ...
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
                error("radia:ngsolve:VectorDeleted", ...
                    "The native NGSolve vector has been deleted.");
            end
        end
    end
end
