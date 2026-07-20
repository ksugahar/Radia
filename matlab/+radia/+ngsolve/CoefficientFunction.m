classdef CoefficientFunction < handle
    %COEFFICIENTFUNCTION Native NGSolve coefficient-function expression.
    %   The numerical expression tree remains in NGSolve. MATLAB stores only
    %   a checked uint64 handle, so arithmetic does not copy field data.

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=protected)
        function obj = CoefficientFunction(nativeHandle)
            obj.NativeHandle = uint64(nativeHandle);
        end
    end

    methods (Static)
        function obj = constant(values)
            arguments
                values double
            end
            obj = radia.ngsolve.CoefficientFunction( ...
                radia.internal.callMex( ...
                'ngsolve.coefficient_function.constant_create', values));
        end

        function obj = fromNativeHandle(nativeHandle)
            obj = radia.ngsolve.CoefficientFunction(nativeHandle);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'ngsolve.coefficient_function.info', obj.NativeHandle);
        end

        function values = evaluate(obj, volPath, points)
            arguments
                obj (1,1) radia.ngsolve.CoefficientFunction
                volPath (1,1) string
                points double
            end
            if ~ismatrix(points) || size(points, 2) < 2 || size(points, 2) > 3
                error("radia:ngsolve:CoefficientFunctionPoints", ...
                    "points must be an N-by-2 or N-by-3 real double matrix.");
            end
            obj.assertAlive();
            values = radia.internal.callMex( ...
                'ngsolve.coefficient_function.evaluate', char(volPath), ...
                obj.NativeHandle, points);
        end

        function result = plus(left, right)
            arguments
                left (1,1) radia.ngsolve.CoefficientFunction
                right (1,1) radia.ngsolve.CoefficientFunction
            end
            result = radia.ngsolve.CoefficientFunction( ...
                radia.internal.callMex('ngsolve.coefficient_function.add', ...
                left.nativeHandle(), right.nativeHandle()));
        end

        function result = minus(left, right)
            arguments
                left (1,1) radia.ngsolve.CoefficientFunction
                right (1,1) radia.ngsolve.CoefficientFunction
            end
            result = radia.ngsolve.CoefficientFunction( ...
                radia.internal.callMex('ngsolve.coefficient_function.subtract', ...
                left.nativeHandle(), right.nativeHandle()));
        end

        function result = mtimes(left, right)
            if isa(left, "radia.ngsolve.CoefficientFunction") && ...
                    isa(right, "radia.ngsolve.CoefficientFunction")
                result = radia.ngsolve.CoefficientFunction( ...
                    radia.internal.callMex('ngsolve.coefficient_function.multiply', ...
                    left.nativeHandle(), right.nativeHandle()));
                return
            end
            if isa(left, "radia.ngsolve.CoefficientFunction")
                coefficient = left;
                scalar = right;
            else
                coefficient = right;
                scalar = left;
            end
            if ~isnumeric(scalar) || ~isscalar(scalar) || ...
                    ~isa(scalar, "double")
                error("radia:ngsolve:CoefficientFunctionScale", ...
                    "A scalar double is required for CoefficientFunction scaling.");
            end
            result = radia.ngsolve.CoefficientFunction( ...
                radia.internal.callMex('ngsolve.coefficient_function.scale', ...
                scalar, coefficient.nativeHandle()));
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'ngsolve.coefficient_function.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:ngsolve:CoefficientFunctionDeleted", ...
                    "The native NGSolve CoefficientFunction has been deleted.");
            end
        end
    end
end
