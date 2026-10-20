classdef HMatrixOperator < handle
    %HMATRIXOPERATOR Experimental resident Gram or B'*G*B application.
    % Numeric leaf schema is an explicit assembly boundary. This class does
    % not construct FE geometry or certify a complete HDiv-MMM solver.
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj=HMatrixOperator(config)
            arguments
                config (1,1) struct
            end
            if isfield(config,'assembly') && strcmp(config.assembly,'hdiv_hacapk_device_dense')
                required={'assembly';'operator';'epsilon';'leaf_size';'eta';'max_rank';'cache_bytes'};
                assert(isequal(sort(fieldnames(config)),sort(required)),'Invalid reusable HDiv configuration');
                obj.Token=radia_cuda_hmatrix('create_hdiv_hacapk_device_dense',config.operator,config.epsilon,config.leaf_size,config.eta,config.max_rank,config.cache_bytes);
                return
            end
            command='create';
            if isfield(config,'assembly')
                if strcmp(config.assembly,'normalized_leaves')
                    command='create_normalized';
                else
                    assert(strcmp(config.assembly,'analytic_constant_block'),'radia:cuda:hmatrix','Unsupported assembly');
                    command='create_analytic_constant_block';
                end
            end
            obj.Token=radia_cuda_hmatrix(command,config);
        end
        function upload(obj,state)
            arguments
                obj (1,1) radia.gpu.HMatrixOperator
                state (:,1) double {mustBeReal,mustBeFinite}
            end
            radia_cuda_hmatrix('upload',obj.Token,state);
        end
        function apply(obj,repetitions)
            arguments
                obj (1,1) radia.gpu.HMatrixOperator
                repetitions (1,1) double {mustBeInteger,mustBeInRange(repetitions,1,10000)} = 1
            end
            % Enqueue all applications; no vector download or host reduction.
            radia_cuda_hmatrix('apply',obj.Token,repetitions);
        end
        function state=download(obj)
            state=radia_cuda_hmatrix('download',obj.Token);
        end
        function info=stats(obj)
            info=radia_cuda_hmatrix('stats',obj.Token);
        end
        function result=solveMaterial(obj,mass,rhs,invChi,tolerance,maxIterations,preconditioner,chebyshevDegree,chebyshevLowerFraction)
            arguments
                obj (1,1) radia.gpu.HMatrixOperator
                mass double {mustBeReal,mustBeFinite}
                rhs (:,1) double {mustBeReal,mustBeFinite}
                invChi (1,1) double {mustBePositive,mustBeFinite}
                tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
                maxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
                preconditioner (1,:) char {mustBeMember(preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial','mass-polynomial-graph','mass-chebyshev','mass-chebyshev-graph'})}='mass-polynomial'
                chebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=8
                chebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}=.1
            end
            % Geometry/charge map remain frozen for this object's lifetime.
            % Mass/RHS are explicit per solve; no hidden mutable source cache.
            assert(chebyshevDegree<=64 && chebyshevLowerFraction<1,'Invalid fixed Chebyshev settings');
            if ismember(preconditioner,{'mass-chebyshev','mass-chebyshev-graph'})
                result=radia_cuda_hmatrix('solve_linear_material',obj.Token,mass,rhs,invChi,tolerance,maxIterations,preconditioner,chebyshevDegree,chebyshevLowerFraction);
            else
                assert(chebyshevDegree==8 && chebyshevLowerFraction==.1,'Chebyshev settings require mass-chebyshev');
                result=radia_cuda_hmatrix('solve_linear_material',obj.Token,mass,rhs,invChi,tolerance,maxIterations,preconditioner);
            end
        end
        function result=solveNonlinearMaterial(obj,material,mass,rhs,initial,options)
            % Shared frozen-tangent Newton and line search; no MATLAB math copy.
            arguments
                obj (1,1) radia.gpu.HMatrixOperator
                material (1,1) radia.gpu.HDivBHMaterial
                mass double {mustBeReal,mustBeFinite}
                rhs (:,1) double {mustBeReal,mustBeFinite}
                initial (:,1) double {mustBeReal,mustBeFinite}
                options (1,1) struct
            end
            assert(issparse(mass),'radia:gpu:HDivNewton','Expected NGSolve sparse reference mass');
            % Native validation owns option bounds, typed handles and epochs.
            % In particular no mass-cholesky -> diagonal substitution occurs.
            result=radia_cuda_hmatrix('solve_nonlinear_material',obj.Token,material.nativeToken(),mass,rhs,initial,options);
        end
        function delete(obj)
            if obj.Token~=0
                radia_cuda_hmatrix('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
end
