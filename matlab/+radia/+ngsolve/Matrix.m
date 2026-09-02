classdef Matrix < handle
    %MATRIX Persistent native NGSolve BaseMatrix handle.
    %   Matrix-vector products stay in the native NGSolve representation.

    properties (SetAccess=private)
        Rows double = 0
        Cols double = 0
        IsSparse logical = false
        IsComplex logical = false
        IsSymmetric logical = false
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = Matrix(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Rows = info.rows;
            obj.Cols = info.cols;
            obj.IsSparse = info.is_sparse;
            obj.IsComplex = info.is_complex;
            obj.IsSymmetric = info.symmetric;
        end
    end

    methods (Static)
        function obj = fromNativeHandle(nativeHandle)
            info = radia.internal.callMex('ngsolve.matrix.info', ...
                nativeHandle);
            obj = radia.ngsolve.Matrix(nativeHandle, info);
        end


        function obj = projected(parent, projection)
            arguments
                parent (1,1) radia.ngsolve.Matrix
                projection (:,:) double
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.projected_create', parent.nativeHandle(), ...
                projection);
            obj = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end

        function obj = reducedBlock(dense, matrices, ranges, scales)
            arguments
                dense (:,:) double
                matrices
                ranges (:,2) double {mustBeInteger, mustBePositive}
                scales (:,1) double
            end
            if size(dense,1) ~= size(dense,2)
                error("radia:ngsolve:ReducedBlockShape", ...
                    "dense must be square.");
            end
            if size(ranges,1) ~= numel(matrices) || ...
                    numel(scales) ~= numel(matrices)
                error("radia:ngsolve:ReducedBlockTerms", ...
                    "matrices, ranges, and scales must have equal lengths.");
            end
            if any(ranges(:,1) > ranges(:,2)) || ...
                    any(ranges(:,2) > size(dense,1))
                error("radia:ngsolve:ReducedBlockRange", ...
                    "ranges must be valid inclusive MATLAB index intervals.");
            end
            handles = zeros(numel(matrices), 1, "uint64");
            for index = 1:numel(matrices)
                if iscell(matrices)
                    matrix = matrices{index};
                else
                    matrix = matrices(index);
                end
                if ~isa(matrix, "radia.ngsolve.Matrix")
                    error("radia:ngsolve:ReducedBlockMatrix", ...
                        "Every embedded term must be a radia.ngsolve.Matrix.");
                end
                handles(index) = matrix.nativeHandle();
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.reduced_block_create', dense, handles, ...
                int32(ranges(:,1)-1), int32(ranges(:,2)), scales(:));
            obj = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.matrix.info', ...
                obj.NativeHandle);
        end

        function triplets = values(obj)
            obj.assertAlive();
            triplets = radia.internal.callMex('ngsolve.matrix.values', ...
                obj.NativeHandle);
        end

        function value = sparse(obj)
            triplets = obj.values();
            value = sparse(triplets.row, triplets.col, triplets.values, ...
                triplets.shape(1), triplets.shape(2));
        end

        function vector = vector(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.vector', obj.NativeHandle);
            vector = radia.ngsolve.Vector.fromNativeHandle(nativeHandle);
        end

        function result = matvec(obj, vector, options)
            arguments
                obj (1,1) radia.ngsolve.Matrix
                vector (1,1) radia.ngsolve.Vector
                options.Transpose (1,1) logical = false
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.matvec', obj.NativeHandle, ...
                vector.nativeHandle(), options.Transpose);
            result = radia.ngsolve.Vector.fromNativeHandle(nativeHandle);
        end

        function matvecInto(obj, input, output, options)
            arguments
                obj (1,1) radia.ngsolve.Matrix
                input (1,1) radia.ngsolve.Vector
                output (1,1) radia.ngsolve.Vector
                options.Transpose (1,1) logical = false
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'ngsolve.matrix.matvec_into', obj.NativeHandle, ...
                input.nativeHandle(), output.nativeHandle(), ...
                options.Transpose);
        end

        function result = inverse(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.inverse', obj.NativeHandle);
            result = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end

        function result = diagonalPreconditioner(obj, relativeFloor)
            arguments
                obj (1,1) radia.ngsolve.Matrix
                relativeFloor (1,1) double {mustBePositive} = 1e-14
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.matrix.diagonal_preconditioner', ...
                obj.NativeHandle, relativeFloor);
            result = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end

        function count = termCount(obj)
            obj.assertAlive();
            count = radia.internal.callMex( ...
                'ngsolve.matrix.term_count', obj.NativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.matrix.destroy', ...
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
                error("radia:ngsolve:MatrixDeleted", ...
                    "The native NGSolve Matrix has been deleted.");
            end
        end
    end
end
