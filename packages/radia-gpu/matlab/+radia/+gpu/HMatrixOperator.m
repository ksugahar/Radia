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
            command='create';
            if isfield(config,'assembly')
                assert(strcmp(config.assembly,'analytic_constant_block'),'radia:cuda:hmatrix','Unsupported assembly');
                command='create_analytic_constant_block';
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
        function delete(obj)
            if obj.Token~=0
                radia_cuda_hmatrix('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
end
