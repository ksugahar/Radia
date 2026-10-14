classdef Solver < handle
    %SOLVER Persistent native NGSolve Krylov solver.
    %   The matrix, optional matrix preconditioner, and iteration state stay
    %   in C++; MATLAB receives only the resulting Vector handle.

    properties (SetAccess=private)
        Method string = ""
        Rows double = 0
        Cols double = 0
        Tolerance double = 0
        MaxSteps double = 0
        IsComplex logical = false
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = Solver(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Method = string(info.method);
            obj.Rows = info.rows;
            obj.Cols = info.cols;
            obj.Tolerance = info.tolerance;
            obj.MaxSteps = info.max_steps;
            obj.IsComplex = info.is_complex;
        end
    end

    methods (Static)
        function obj = create(matrix, method, options)
            arguments
                matrix (1,1) radia.ngsolve.Matrix
                method (1,1) string
                options.Tolerance (1,1) double = 1e-8
                options.MaxSteps (1,1) double {mustBeInteger, mustBePositive} = 1000
                options.Preconditioner = []
            end
            if isempty(options.Preconditioner)
                nativeHandle = radia.internal.callMex( ...
                    'ngsolve.solver.create', matrix.nativeHandle(), ...
                    char(method), options.Tolerance, options.MaxSteps);
            else
                if ~isa(options.Preconditioner, "radia.ngsolve.Matrix")
                    error("radia:ngsolve:SolverPreconditioner", ...
                        "Preconditioner must be a radia.ngsolve.Matrix.");
                end
                nativeHandle = radia.internal.callMex( ...
                    'ngsolve.solver.create', matrix.nativeHandle(), ...
                    char(method), options.Tolerance, options.MaxSteps, ...
                    options.Preconditioner.nativeHandle());
            end
            info = radia.internal.callMex('ngsolve.solver.info', nativeHandle);
            obj = radia.ngsolve.Solver(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.solver.info', ...
                obj.NativeHandle);
        end

        function solution = solve(obj, rhs)
            arguments
                obj (1,1) radia.ngsolve.Solver
                rhs (1,1) radia.ngsolve.Vector
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.solver.solve', obj.NativeHandle, rhs.nativeHandle());
            solution = radia.ngsolve.Vector.fromNativeHandle(nativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.solver.destroy', ...
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
                error("radia:ngsolve:SolverDeleted", ...
                    "The native NGSolve Solver has been deleted.");
            end
        end
    end
end
