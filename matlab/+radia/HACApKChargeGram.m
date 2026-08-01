classdef HACApKChargeGram < handle
    %HACAPKCHARGEGRAM Stateful HDiv charge Gram H-matrix.

    properties (SetAccess = private)
        NDOF = 0
    end

    properties (Access = private)
        NativeHandle = uint64(0)
    end

    methods
        function obj = HACApKChargeGram(centroids, measures, selfEnergy)
            if nargin == 0
                return
            end
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_monopole', double(centroids), ...
                double(measures), double(selfEnergy));
            obj.refreshNDOF();
        end

        function ok = build(obj, options)
            arguments
                obj
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
            end
            obj.assertAlive();
            ok = radia.internal.callMex( ...
                'hacapk.charge_gram.build', obj.NativeHandle, ...
                options.AcaEps, options.LeafSize, options.Eta, ...
                options.MaxRank, options.PrintLevel);
            if ok
                obj.NDOF = radia.internal.callMex( ...
                    'hacapk.charge_gram.info', obj.NativeHandle).n_dof;
            end
        end

        function y = matvec(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec', x);
        end

        function y = matvecTranspose(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec_transpose', x);
        end

        function y = matvecSym(obj, x)
            y = obj.apply('hacapk.charge_gram.matvec_sym', x);
        end

        function value = entry(obj, i, j)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'hacapk.charge_gram.entry', obj.NativeHandle, i, j);
        end

        function block = hexVolumeSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (27,3) double
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_volume_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function block = hexFaceSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (9,3) double
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_face_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function matrix = hexChargeGramDirectionalDerivative( ...
                obj, cellNodeVelocity, faceNodeVelocity)
            if ndims(cellNodeVelocity) ~= 3 || size(cellNodeVelocity,2) ~= 27 || ...
                    size(cellNodeVelocity,3) ~= 3
                error("radia:HACApKChargeGram:HexCellVelocityShape", ...
                    "cellNodeVelocity must have shape nCell-by-27-by-3.");
            end
            if ndims(faceNodeVelocity) ~= 3 || size(faceNodeVelocity,2) ~= 9 || ...
                    size(faceNodeVelocity,3) ~= 3
                error("radia:HACApKChargeGram:HexFaceVelocityShape", ...
                    "faceNodeVelocity must have shape nFace-by-9-by-3.");
            end
            obj.assertAlive();
            matrix = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_directional_derivative', ...
                obj.NativeHandle, double(cellNodeVelocity), ...
                double(faceNodeVelocity));
        end

        function block = tetVolumeSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (4,3) double
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.tet_volume_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function block = tetFaceSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (3,3) double
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.tet_face_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function matrix = tetChargeGramDirectionalDerivative( ...
                obj, cellVertexVelocity, faceVertexVelocity)
            if ndims(cellVertexVelocity) ~= 3 || size(cellVertexVelocity,2) ~= 4 || ...
                    size(cellVertexVelocity,3) ~= 3
                error("radia:HACApKChargeGram:TetCellVelocityShape", ...
                    "cellVertexVelocity must have shape nCell-by-4-by-3.");
            end
            if ndims(faceVertexVelocity) ~= 3 || size(faceVertexVelocity,2) ~= 3 || ...
                    size(faceVertexVelocity,3) ~= 3
                error("radia:HACApKChargeGram:TetFaceVelocityShape", ...
                    "faceVertexVelocity must have shape nFace-by-3-by-3.");
            end
            obj.assertAlive();
            matrix = radia.internal.callMex( ...
                'hacapk.charge_gram.tet_directional_derivative', ...
                obj.NativeHandle, double(cellVertexVelocity), ...
                double(faceVertexVelocity));
        end

        function rates = tetChargeMapRowDirectionalRates( ...
                obj, cellVertexVelocity, faceVertexVelocity)
            if ndims(cellVertexVelocity) ~= 3 || size(cellVertexVelocity,2) ~= 4 || ...
                    size(cellVertexVelocity,3) ~= 3
                error("radia:HACApKChargeGram:TetCellVelocityShape", ...
                    "cellVertexVelocity must have shape nCell-by-4-by-3.");
            end
            if ndims(faceVertexVelocity) ~= 3 || size(faceVertexVelocity,2) ~= 3 || ...
                    size(faceVertexVelocity,3) ~= 3
                error("radia:HACApKChargeGram:TetFaceVelocityShape", ...
                    "faceVertexVelocity must have shape nFace-by-3-by-3.");
            end
            obj.assertAlive();
            rates = radia.internal.callMex( ...
                'hacapk.charge_gram.tet_charge_map_row_directional_rates', ...
                obj.NativeHandle, double(cellVertexVelocity), ...
                double(faceVertexVelocity));
        end

        function block = wedgeVolumeSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (18,3) double
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.wedge_volume_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function block = wedgeFaceSelfBlockDirectionalDerivative( ...
                obj, host, nodeVelocity)
            arguments
                obj
                host (1,1) double
                nodeVelocity (:,3) double
            end
            if ~ismember(size(nodeVelocity,1), [6, 9])
                error("radia:HACApKChargeGram:WedgeFaceShape", ...
                    "nodeVelocity must have 6 or 9 rows.");
            end
            obj.assertAlive();
            block = radia.internal.callMex( ...
                'hacapk.charge_gram.wedge_face_self_block_directional_derivative', ...
                obj.NativeHandle, host, nodeVelocity);
        end

        function matrix = wedgeChargeGramDirectionalDerivative( ...
                obj, cellNodeVelocity, faceNodeVelocity)
            if ndims(cellNodeVelocity) ~= 3 || size(cellNodeVelocity,2) ~= 18 || ...
                    size(cellNodeVelocity,3) ~= 3
                error("radia:HACApKChargeGram:WedgeCellVelocityShape", ...
                    "cellNodeVelocity must have shape nCell-by-18-by-3.");
            end
            if ndims(faceNodeVelocity) ~= 3 || size(faceNodeVelocity,2) ~= 9 || ...
                    size(faceNodeVelocity,3) ~= 3
                error("radia:HACApKChargeGram:WedgeFaceVelocityShape", ...
                    "faceNodeVelocity must have shape nFace-by-9-by-3.");
            end
            obj.assertAlive();
            matrix = radia.internal.callMex( ...
                'hacapk.charge_gram.wedge_directional_derivative', ...
                obj.NativeHandle, double(cellNodeVelocity), ...
                double(faceNodeVelocity));
        end

        function derivative = directionalDerivativeOperator( ...
                obj, family, cellVelocity, faceVelocity, options)
            arguments
                obj
                family (1,1) string
                cellVelocity double
                faceVelocity double
                options.AcaEps (1,1) double {mustBePositive} = 1e-8
                options.LeafSize (1,1) double {mustBeInteger,mustBePositive} = 32
                options.Eta (1,1) double {mustBePositive} = 2
            end
            family = lower(family);
            if ~ismember(family, ["hex", "tet", "wedge"])
                error("radia:HACApKChargeGram:DerivativeFamily", ...
                    "family must be hex, tet, or wedge.");
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.directional_derivative_operator', ...
                obj.NativeHandle, char(family), double(cellVelocity), ...
                double(faceVelocity), options.AcaEps, options.LeafSize, ...
                options.Eta);
            derivative = ...
                radia.HACApKChargeGramDerivative.fromNativeHandle(nativeHandle);
        end

        function values = directionalDerivativeContractions( ...
                obj, family, cellVelocity, faceVelocity, left, right)
            arguments
                obj
                family (1,1) string
                cellVelocity double
                faceVelocity double
                left (:,1) double
                right (:,1) double
            end
            family = lower(family);
            if ~ismember(family, ["hex", "tet", "wedge"])
                error("radia:HACApKChargeGram:DerivativeFamily", ...
                    "family must be hex, tet, or wedge.");
            end
            if ndims(cellVelocity) ~= 4 || ndims(faceVelocity) ~= 4 || ...
                    size(cellVelocity,1) ~= size(faceVelocity,1) || ...
                    size(cellVelocity,4) ~= 3 || size(faceVelocity,4) ~= 3
                error("radia:HACApKChargeGram:DerivativeBatchShape", ...
                    ["Velocity arrays must have shape " ...
                     "nMode-by-nHost-by-nNode-by-3 with equal mode counts."]);
            end
            expected = struct("hex",[27,9],"tet",[4,3],"wedge",[18,9]);
            nodeCounts = expected.(family);
            if size(cellVelocity,3) ~= nodeCounts(1) || ...
                    size(faceVelocity,3) ~= nodeCounts(2)
                error("radia:HACApKChargeGram:DerivativeBatchShape", ...
                    "Velocity node counts do not match the selected family.");
            end
            obj.assertAlive();
            values = radia.internal.callMex( ...
                'hacapk.charge_gram.directional_derivative_contractions', ...
                obj.NativeHandle, char(family), double(cellVelocity), ...
                double(faceVelocity), double(left), double(right));
        end

        function values = directionalDerivativeContractionsMany( ...
                obj, family, cellVelocity, faceVelocity, left, right)
            arguments
                obj
                family (1,1) string
                cellVelocity double
                faceVelocity double
                left (:,:) double
                right (:,1) double
            end
            family = lower(family);
            if ~ismember(family, ["hex", "tet", "wedge"])
                error("radia:HACApKChargeGram:DerivativeFamily", ...
                    "family must be hex, tet, or wedge.");
            end
            if ndims(cellVelocity) ~= 4 || ndims(faceVelocity) ~= 4 || ...
                    size(cellVelocity,1) ~= size(faceVelocity,1) || ...
                    size(cellVelocity,4) ~= 3 || size(faceVelocity,4) ~= 3
                error("radia:HACApKChargeGram:DerivativeBatchShape", ...
                    ["Velocity arrays must have shape " ...
                     "nMode-by-nHost-by-nNode-by-3 with equal mode counts."]);
            end
            expected = struct("hex",[27,9],"tet",[4,3],"wedge",[18,9]);
            nodeCounts = expected.(family);
            if size(cellVelocity,3) ~= nodeCounts(1) || ...
                    size(faceVelocity,3) ~= nodeCounts(2)
                error("radia:HACApKChargeGram:DerivativeBatchShape", ...
                    "Velocity node counts do not match the selected family.");
            end
            obj.assertAlive();
            values = radia.internal.callMex( ...
                'hacapk.charge_gram.directional_derivative_contractions_many', ...
                obj.NativeHandle, char(family), double(cellVelocity), ...
                double(faceVelocity), double(left), double(right));
        end

        function info = info(obj)
            obj.assertAlive();
            info = radia.internal.callMex( ...
                'hacapk.charge_gram.info', obj.NativeHandle);
        end

        function result = hexStateCheck(obj)
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_state_check', obj.NativeHandle);
        end

        function result = hexStoredNodes(obj)
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_stored_nodes', obj.NativeHandle);
        end

        function result = hexStateBreakdown(obj)
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.hex_state_breakdown', obj.NativeHandle);
        end

        function configureChargeMap(obj, indptr, indices, data, nFace)
            obj.assertAlive();
            radia.internal.callMex('hacapk.charge_gram.configure_charge_map', ...
                obj.NativeHandle, int32(indptr), int32(indices), double(data), nFace);
        end

        function configureVectorChargeMap(obj, indptr, indices, data, nFace, nComponents)
            arguments
                obj
                indptr
                indices
                data
                nFace (1,1) double
                nComponents (1,1) double = 3
            end
            obj.assertAlive();
            radia.internal.callMex('hacapk.charge_gram.configure_vector_charge_map', ...
                obj.NativeHandle, int32(indptr), int32(indices), double(data), nFace, nComponents);
        end

        function configureMassMatrix(obj, rows, cols, values, nFace)
            obj.assertAlive();
            radia.internal.callMex('hacapk.charge_gram.configure_mass_matrix', ...
                obj.NativeHandle, int32(rows), int32(cols), double(values), nFace);
        end

        function configureGeometryMassMatrix(obj, rows, cols, values, nFace)
            obj.assertAlive();
            radia.internal.callMex('hacapk.charge_gram.configure_geometry_mass_matrix', ...
                obj.NativeHandle, int32(rows), int32(cols), double(values), nFace);
        end

        function configureMassMatrixNGSolve(obj, matrix)
            arguments
                obj (1,1) radia.HACApKChargeGram
                matrix (1,1) radia.ngsolve.Matrix
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'hacapk.charge_gram.configure_mass_matrix_ngsolve', ...
                obj.NativeHandle, matrix.nativeHandle());
        end

        function configureGeometryMassMatrixNGSolve(obj, matrix)
            arguments
                obj (1,1) radia.HACApKChargeGram
                matrix (1,1) radia.ngsolve.Matrix
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'hacapk.charge_gram.configure_geometry_mass_matrix_ngsolve', ...
                obj.NativeHandle, matrix.nativeHandle());
        end

        function changed = restoreGeometryMassMatrix(obj)
            obj.assertAlive();
            changed = radia.internal.callMex( ...
                'hacapk.charge_gram.restore_geometry_mass_matrix', obj.NativeHandle);
        end

        function setConfiguredConstraints(obj, dofs, options)
            arguments
                obj
                dofs
                options.PreserveExisting (1,1) logical = false
            end
            obj.assertAlive();
            radia.internal.callMex( ...
                'hacapk.charge_gram.set_configured_constraints', ...
                obj.NativeHandle, int32(dofs), options.PreserveExisting);
        end

        function result = operatorInfo(obj)
            obj.assertAlive();
            result = radia.internal.callMex('hacapk.charge_gram.operator_info', obj.NativeHandle);
        end

        function matrix = demagMatrix(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.demag_matrix', obj.NativeHandle);
            matrix = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end

        function y = applyConfiguredDemag(obj, x, options)
            arguments
                obj
                x
                options.Symmetric (1,1) logical = true
            end
            obj.assertAlive();
            y = radia.internal.callMex('hacapk.charge_gram.demag_apply', ...
                obj.NativeHandle, double(x), options.Symmetric);
        end

        function y = applyConfiguredGeometryMass(obj, x)
            obj.assertAlive();
            y = radia.internal.callMex('hacapk.charge_gram.geometry_mass_apply', ...
                obj.NativeHandle, double(x));
        end

        function y = applyConfiguredMassRiesz(obj, rhs)
            obj.assertAlive();
            y = radia.internal.callMex('hacapk.charge_gram.mass_riesz', ...
                obj.NativeHandle, double(rhs));
        end

        function y = applyConfiguredLinearMaterialOperator( ...
                obj, invChi, x, options)
            arguments
                obj
                invChi (1,1) double
                x
                options.RespectConstraints (1,1) logical = true
            end
            obj.assertAlive();
            y = radia.internal.callMex( ...
                'hacapk.charge_gram.apply_configured_linear_material_operator', ...
                obj.NativeHandle, invChi, double(x), ...
                options.RespectConstraints);
        end

        function y = applyConfiguredLinearMaterialOperatorMany( ...
                obj, invChi, x, options)
            arguments
                obj
                invChi (1,1) double
                x (:,:) double
                options.RespectConstraints (1,1) logical = true
            end
            obj.assertAlive();
            y = radia.internal.callMex( ...
                ['hacapk.charge_gram.' ...
                 'apply_configured_linear_material_operator_many'], ...
                obj.NativeHandle, invChi, double(x), ...
                options.RespectConstraints);
        end

        function result = reduceConfiguredCandidateSchur( ...
                obj, invChi, candidateDofs, rhs, state, responseMatrix, ...
                adjoints, options)
            arguments
                obj
                invChi (1,1) double
                candidateDofs
                rhs
                state
                responseMatrix (:,:) double
                adjoints (:,:) double
                options.Tol (1,1) double = 1e-9
                options.MaxIt (1,1) double = 5000
                options.SolveBatchSize (1,1) double = 64
                options.MassRiesz (1,1) logical = true
            end
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.reduce_configured_candidate_schur', ...
                obj.NativeHandle, invChi, int32(candidateDofs), double(rhs), ...
                double(state), double(responseMatrix), double(adjoints), ...
                options.Tol, options.MaxIt, options.SolveBatchSize, ...
                options.MassRiesz);
        end

        function result = solveConfiguredLinearMaterial(obj, invChi, rhs, options)
            arguments
                obj
                invChi (1,1) double
                rhs
                options.Tol (1,1) double = 1e-8
                options.MaxIt (1,1) double = 5000
                options.Symmetric (1,1) logical = true
                options.X0 = []
            end
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.solve_configured_linear_material', ...
                obj.NativeHandle, invChi, double(rhs), options.Tol, options.MaxIt, ...
                options.Symmetric, double(options.X0));
        end

        function result = solveConfiguredLinearMaterialAutoPrec(obj, invChi, rhs, options)
            arguments
                obj
                invChi (1,1) double
                rhs
                options.Tol (1,1) double = 1e-9
                options.MaxIt (1,1) double = 5000
                options.X0 = []
            end
            obj.assertAlive();
            result = radia.internal.callMex( ...
                'hacapk.charge_gram.solve_configured_linear_material_auto_prec', ...
                obj.NativeHandle, invChi, double(rhs), options.Tol, options.MaxIt, ...
                true, double(options.X0));
        end

        function result = solveConfiguredLinearMaterialAutoPrecMany( ...
                obj, invChi, rhs, options)
            arguments
                obj
                invChi (1,1) double
                rhs (:,:) double
                options.Tol (1,1) double = 1e-9
                options.MaxIt (1,1) double = 5000
                options.ClusterCoarseSize (1,1) double = 0
                options.ClusterDeflationSize (1,1) double = 0
                options.RecycleSize (1,1) double = 0
                options.MassRiesz (1,1) logical = false
                options.X0 = []
            end
            obj.assertAlive();
            result = radia.internal.callMex( ...
                ['hacapk.charge_gram.' ...
                 'solve_configured_linear_material_auto_prec_many'], ...
                obj.NativeHandle, invChi, double(rhs), options.Tol, ...
                options.MaxIt, options.ClusterCoarseSize, ...
                options.ClusterDeflationSize, options.RecycleSize, ...
                options.MassRiesz, double(options.X0));
        end

        function rows = configuredFieldFunctionalRows( ...
                obj, observations, weights)
            arguments
                obj
                observations (:,3) double
                weights double
            end
            obj.assertAlive();
            rows = radia.internal.callMex( ...
                'hacapk.charge_gram.configured_field_functional_rows', ...
                obj.NativeHandle, double(observations), double(weights));
        end

        function rows = configuredFieldFunctionalRowsDirectionalDerivative( ...
                obj, observations, weights, cellVertexVelocity, ...
                faceVertexVelocity)
            arguments
                obj
                observations (:,3) double
                weights double
                cellVertexVelocity double
                faceVertexVelocity double
            end
            obj.assertAlive();
            rows = radia.internal.callMex( ...
                ['hacapk.charge_gram.' ...
                 'configured_field_functional_rows_directional_derivative'], ...
                obj.NativeHandle, double(observations), double(weights), ...
                double(cellVertexVelocity), double(faceVertexVelocity));
        end

        function evaluator = createFieldEvaluator(obj, magnetization, options)
            arguments
                obj
                magnetization
                options.LeafSize (1,1) double = 32
                options.Theta (1,1) double = 0.05
                options.TreeMinSources (1,1) double = 256
                options.AutoMinWork (1,1) double = 500000000
                options.TreeRelativeTolerance (1,1) double = 1e-5
                options.ProbeCount (1,1) double = 16
            end
            obj.assertAlive();
            evaluator = radia.HDivFieldEvaluator.fromChargeGram(obj, magnetization, ...
                LeafSize=options.LeafSize, Theta=options.Theta, ...
                TreeMinSources=options.TreeMinSources, AutoMinWork=options.AutoMinWork, ...
                TreeRelativeTolerance=options.TreeRelativeTolerance, ProbeCount=options.ProbeCount);
        end

        function evaluator = createPlanarFieldEvaluator(obj, magnetization)
            obj.assertAlive();
            evaluator = radia.PlanarFieldEvaluator.fromChargeGram(obj, magnetization);
        end

        function result = stats(obj)
            obj.assertAlive();
            result = radia.internal.callMex('hacapk.charge_gram.stats', obj.NativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'hacapk.charge_gram.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Static)
        function obj = from_sampled_laplace(points, weights, kernelEpsilon, options)
            arguments
                points double
                weights double
                kernelEpsilon (1,1) double
                options.AcaEps (1,1) double = 1e-11
                options.LeafSize (1,1) double = 64
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_sampled_laplace', ...
                double(points), double(weights), kernelEpsilon);
            obj.finishConstruction(options);
        end

        function obj = from_sampled_planar_log( ...
                points, weights, kernelEpsilon, referenceLength, options)
            arguments
                points double
                weights double
                kernelEpsilon (1,1) double
                referenceLength (1,1) double = 1
                options.AcaEps (1,1) double = 1e-11
                options.LeafSize (1,1) double = 64
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_sampled_planar_log', ...
                double(points), double(weights), kernelEpsilon, referenceLength);
            obj.finishConstruction(options);
        end

        function obj = from_local_polynomials( ...
                cellVerts, nElements, chargeHost, polynomialCoefficients, ...
                polynomialExponents, referenceTetPoints, referenceTetWeights, options)
            arguments
                cellVerts double
                nElements (1,1) double
                chargeHost
                polynomialCoefficients double
                polynomialExponents double
                referenceTetPoints double
                referenceTetWeights double
                options.AcaEps (1,1) double = 1e-6
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_local_polynomials', ...
                double(cellVerts), nElements, int32(chargeHost), ...
                double(polynomialCoefficients), double(polynomialExponents), ...
                double(referenceTetPoints), double(referenceTetWeights));
            obj.finishConstruction(options);
        end

        function obj = from_high_order_tet( ...
                cellVerts, faceVerts, nElements, chargeHost, chargeKind, ...
                chargeExponents, referenceTetPoints, referenceTetWeights, ...
                referenceTrianglePoints, referenceTriangleWeights, options)
            arguments
                cellVerts double
                faceVerts double
                nElements (1,1) double
                chargeHost
                chargeKind
                chargeExponents
                referenceTetPoints double
                referenceTetWeights double
                referenceTrianglePoints double
                referenceTriangleWeights double
                options.LowTetPoints double = zeros(0,3)
                options.LowTetWeights double = zeros(0,1)
                options.LowTrianglePoints double = zeros(0,2)
                options.LowTriangleWeights double = zeros(0,1)
                options.FarFactor (1,1) double = 1e30
                options.InnerTetPoints double = zeros(0,3)
                options.InnerTetWeights double = zeros(0,1)
                options.InnerTrianglePoints double = zeros(0,2)
                options.InnerTriangleWeights double = zeros(0,1)
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_high_order_tet', ...
                double(cellVerts), double(faceVerts), nElements, ...
                int32(chargeHost), int32(chargeKind), int32(chargeExponents), ...
                double(referenceTetPoints), double(referenceTetWeights), ...
                double(referenceTrianglePoints), double(referenceTriangleWeights), ...
                double(options.LowTetPoints), double(options.LowTetWeights), ...
                double(options.LowTrianglePoints), double(options.LowTriangleWeights), ...
                options.FarFactor, double(options.InnerTetPoints), ...
                double(options.InnerTetWeights), double(options.InnerTrianglePoints), ...
                double(options.InnerTriangleWeights), int32(options.ImageMasks), ...
                double(options.ImageSigns));
            obj.finishConstruction(options);
        end

        function obj = from_curved_high_order_tet( ...
                cellNodes, faceNodes, cellVertices, faceVertices, nElements, ...
                chargeHost, chargeKind, chargeExponents, referenceTetPoints, ...
                referenceTetWeights, referenceTrianglePoints, ...
                referenceTriangleWeights, curvePoints, curveWeights, options)
            arguments
                cellNodes double
                faceNodes double
                cellVertices
                faceVertices
                nElements (1,1) double
                chargeHost
                chargeKind
                chargeExponents
                referenceTetPoints double
                referenceTetWeights double
                referenceTrianglePoints double
                referenceTriangleWeights double
                curvePoints double
                curveWeights double
                options.CurveOrder (1,1) double = 2
                options.LowTetPoints double = zeros(0,3)
                options.LowTetWeights double = zeros(0,1)
                options.LowTrianglePoints double = zeros(0,2)
                options.LowTriangleWeights double = zeros(0,1)
                options.FarFactor (1,1) double = 1e30
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.ReferenceDensity (1,1) logical = false
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_curved_high_order_tet', ...
                double(cellNodes), double(faceNodes), int32(cellVertices), ...
                int32(faceVertices), nElements, options.CurveOrder, ...
                int32(chargeHost), int32(chargeKind), int32(chargeExponents), ...
                double(referenceTetPoints), double(referenceTetWeights), ...
                double(referenceTrianglePoints), double(referenceTriangleWeights), ...
                double(curvePoints), double(curveWeights), ...
                double(options.LowTetPoints), double(options.LowTetWeights), ...
                double(options.LowTrianglePoints), double(options.LowTriangleWeights), ...
                options.FarFactor, int32(options.ImageMasks), ...
                double(options.ImageSigns), options.ReferenceDensity);
            obj.finishConstruction(options);
        end

        function obj = from_hex( ...
                hexCellNodes, quadFaceNodes, nElements, nBoundaryFaces, ...
                chargeHost, chargeKind, chargeExponents, symmetricTetPoints, ...
                symmetricTetWeights, symmetricTrianglePoints, ...
                symmetricTriangleWeights, outerPoints, outerWeights, ...
                innerPoints, innerWeights, farTetPoints, farTetWeights, ...
                farTrianglePoints, farTriangleWeights, options)
            arguments
                hexCellNodes double
                quadFaceNodes double
                nElements (1,1) double
                nBoundaryFaces (1,1) double
                chargeHost
                chargeKind
                chargeExponents
                symmetricTetPoints double
                symmetricTetWeights double
                symmetricTrianglePoints double
                symmetricTriangleWeights double
                outerPoints double
                outerWeights double
                innerPoints double
                innerWeights double
                farTetPoints double
                farTetWeights double
                farTrianglePoints double
                farTriangleWeights double
                options.NearGrade (1,1) double = 1.5
                options.FarInnerFactor (1,1) double = 1.5
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_hex', double(hexCellNodes), ...
                double(quadFaceNodes), nElements, nBoundaryFaces, ...
                int32(chargeHost), int32(chargeKind), int32(chargeExponents), ...
                double(symmetricTetPoints), double(symmetricTetWeights), ...
                double(symmetricTrianglePoints), double(symmetricTriangleWeights), ...
                double(outerPoints), double(outerWeights), ...
                double(innerPoints), double(innerWeights), ...
                double(farTetPoints), double(farTetWeights), ...
                double(farTrianglePoints), double(farTriangleWeights), ...
                options.NearGrade, options.FarInnerFactor, ...
                int32(options.ImageMasks), double(options.ImageSigns));
            obj.finishConstruction(options);
        end

        function obj = from_wedge( ...
                wedgeCellNodes, faceNodes, faceType, nElements, nBoundaryFaces, ...
                chargeHost, chargeKind, chargeExponents, symmetricTetPoints, ...
                symmetricTetWeights, symmetricTrianglePoints, ...
                symmetricTriangleWeights, fieldTrianglePoints, ...
                fieldTriangleWeights, outerPoints, outerWeights, innerPoints, ...
                innerWeights, farTetPoints, farTetWeights, farTrianglePoints, ...
                farTriangleWeights, options)
            arguments
                wedgeCellNodes double
                faceNodes double
                faceType
                nElements (1,1) double
                nBoundaryFaces (1,1) double
                chargeHost
                chargeKind
                chargeExponents
                symmetricTetPoints double
                symmetricTetWeights double
                symmetricTrianglePoints double
                symmetricTriangleWeights double
                fieldTrianglePoints double
                fieldTriangleWeights double
                outerPoints double
                outerWeights double
                innerPoints double
                innerWeights double
                farTetPoints double
                farTetWeights double
                farTrianglePoints double
                farTriangleWeights double
                options.NearGrade (1,1) double = 0.6
                options.FarInnerFactor (1,1) double = 1.5
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.AcaEps (1,1) double = 1e-12
                options.LeafSize (1,1) double = 64
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_wedge', double(wedgeCellNodes), ...
                double(faceNodes), int32(faceType), nElements, nBoundaryFaces, ...
                int32(chargeHost), int32(chargeKind), int32(chargeExponents), ...
                double(symmetricTetPoints), double(symmetricTetWeights), ...
                double(symmetricTrianglePoints), double(symmetricTriangleWeights), ...
                double(fieldTrianglePoints), double(fieldTriangleWeights), ...
                double(outerPoints), double(outerWeights), ...
                double(innerPoints), double(innerWeights), ...
                double(farTetPoints), double(farTetWeights), ...
                double(farTrianglePoints), double(farTriangleWeights), ...
                options.NearGrade, options.FarInnerFactor, ...
                int32(options.ImageMasks), double(options.ImageSigns));
            obj.finishConstruction(options);
        end

        function obj = from_planar_2d( ...
                dim2, geometryOrder, cellMap, cellType, edgeMap, nElements, ...
                nBoundaryEdges, chargeHost, chargeKind, chargeExponents, ...
                symmetricTrianglePoints, symmetricTriangleWeights, ...
                quadrilateralPoints, quadrilateralWeights, edgePoints, ...
                edgeWeights, innerPoints, innerWeights, farTrianglePoints, ...
                farTriangleWeights, options)
            arguments
                dim2 (1,1) double
                geometryOrder (1,1) double
                cellMap double
                cellType
                edgeMap double
                nElements (1,1) double
                nBoundaryEdges (1,1) double
                chargeHost
                chargeKind
                chargeExponents
                symmetricTrianglePoints double
                symmetricTriangleWeights double
                quadrilateralPoints double
                quadrilateralWeights double
                edgePoints double
                edgeWeights double
                innerPoints double
                innerWeights double
                farTrianglePoints double
                farTriangleWeights double
                options.NearGrade (1,1) double = 0.6
                options.FarInnerFactor (1,1) double = 1.5
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.AcaEps (1,1) double = 1e-12
                options.LeafSize (1,1) double = 64
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_planar_2d', dim2, geometryOrder, ...
                double(cellMap), int32(cellType), double(edgeMap), ...
                nElements, nBoundaryEdges, int32(chargeHost), ...
                int32(chargeKind), int32(chargeExponents), ...
                double(symmetricTrianglePoints), ...
                double(symmetricTriangleWeights), ...
                double(quadrilateralPoints), double(quadrilateralWeights), ...
                double(edgePoints), double(edgeWeights), ...
                double(innerPoints), double(innerWeights), ...
                double(farTrianglePoints), double(farTriangleWeights), ...
                options.NearGrade, options.FarInnerFactor, ...
                int32(options.ImageMasks), double(options.ImageSigns));
            obj.finishConstruction(options);
        end

        function obj = from_curved_polytope( ...
                cellCurvedNodes, cellSubtetOffsets, cellCentroids, ...
                cellMeasures, faceCurvedNodes, faceSubtriangleOffsets, ...
                faceCentroids, faceMeasures, referenceTetPoints, ...
                referenceTetWeights, referenceTrianglePoints, ...
                referenceTriangleWeights, curvePoints, curveWeights, ...
                nElements, options)
            arguments
                cellCurvedNodes double
                cellSubtetOffsets
                cellCentroids double
                cellMeasures double
                faceCurvedNodes double
                faceSubtriangleOffsets
                faceCentroids double
                faceMeasures double
                referenceTetPoints double
                referenceTetWeights double
                referenceTrianglePoints double
                referenceTriangleWeights double
                curvePoints double
                curveWeights double
                nElements (1,1) double
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_curved_polytope', ...
                double(cellCurvedNodes), int32(cellSubtetOffsets), ...
                double(cellCentroids), double(cellMeasures), ...
                double(faceCurvedNodes), int32(faceSubtriangleOffsets), ...
                double(faceCentroids), double(faceMeasures), ...
                double(referenceTetPoints), double(referenceTetWeights), ...
                double(referenceTrianglePoints), ...
                double(referenceTriangleWeights), double(curvePoints), ...
                double(curveWeights), nElements);
            obj.finishConstruction(options);
        end

        function obj = from_analytic_tet(cellVerts, faceVerts, nElements, options)
            arguments
                cellVerts double
                faceVerts double
                nElements (1,1) double
                options.NearFactor (1,1) double = 1e30
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.FarQuadrature (1,1) double = 0
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_analytic_tet', ...
                double(cellVerts), double(faceVerts), nElements, ...
                options.NearFactor, int32(options.ImageMasks), ...
                double(options.ImageSigns), options.FarQuadrature);
            obj.finishConstruction(options);
        end

        function obj = from_analytic_polytope( ...
                cellTriangles, cellTriangleOffsets, cellCentroids, cellMeasures, ...
                faceTriangles, faceTriangleOffsets, faceCentroids, faceMeasures, ...
                nElements, options)
            arguments
                cellTriangles double
                cellTriangleOffsets
                cellCentroids double
                cellMeasures double
                faceTriangles double
                faceTriangleOffsets
                faceCentroids double
                faceMeasures double
                nElements (1,1) double
                options.NearFactor (1,1) double = 1e30
                options.ImageMasks = int32.empty
                options.ImageSigns double = zeros(0,1)
                options.FarQuadrature (1,1) double = 0
                options.AcaEps (1,1) double = 1e-4
                options.LeafSize (1,1) double = 32
                options.Eta (1,1) double = 2
                options.MaxRank (1,1) double = 200
                options.PrintLevel (1,1) double = 0
                options.Build (1,1) logical = true
            end
            obj = radia.HACApKChargeGram();
            obj.NativeHandle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_analytic_polytope', ...
                double(cellTriangles), int32(cellTriangleOffsets), ...
                double(cellCentroids), double(cellMeasures), ...
                double(faceTriangles), int32(faceTriangleOffsets), ...
                double(faceCentroids), double(faceMeasures), nElements, ...
                options.NearFactor, int32(options.ImageMasks), ...
                double(options.ImageSigns), options.FarQuadrature);
            obj.finishConstruction(options);
        end
    end

    methods (Access = private)
        function finishConstruction(obj, options)
            obj.refreshNDOF();
            if options.Build
                ok = obj.build(AcaEps=options.AcaEps, ...
                    LeafSize=options.LeafSize, Eta=options.Eta, ...
                    MaxRank=options.MaxRank, PrintLevel=options.PrintLevel);
                if ~ok
                    delete(obj);
                    error("radia:HACApKChargeGram:BuildFailed", ...
                        "The native HACApK charge-Gram build failed.");
                end
            end
        end

        function refreshNDOF(obj)
            obj.NDOF = radia.internal.callMex( ...
                'hacapk.charge_gram.info', obj.NativeHandle).n_dof;
        end

        function y = apply(obj, command, x)
            obj.assertAlive();
            y = radia.internal.callMex(command, obj.NativeHandle, double(x));
        end

        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HACApKChargeGram:Deleted", ...
                    "The native HACApK charge-Gram manager has been deleted.");
            end
        end
    end
end
