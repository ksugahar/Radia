classdef HDivBHMaterial < handle
    %HDIVBHMATERIAL Persistent shared CUDA inverse-BH state.
    % E and weights come from NGSolve material quadrature in its numbering.
    % This owner does not assemble FE geometry or select a CPU fallback.
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj=HDivBHMaterial(E,weights,law)
            arguments
                E double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite,mustBePositive}
                law (1,1) struct
            end
            assert(issparse(E),'radia:gpu:BHMaterial','Expected NGSolve sparse evaluation matrix');
            obj.Token=radia_cuda_hmatrix('create_bh_material',E,weights,law);
        end
        function result=evaluate(obj,state)
            arguments
                obj (1,1) radia.gpu.HDivBHMaterial
                state (:,1) double {mustBeReal,mustBeFinite}
            end
            result=radia_cuda_hmatrix('evaluate_bh_material_context',obj.nativeToken(),state,[]);
        end
        function info=stats(obj)
            info=radia_cuda_hmatrix('bh_material_stats',obj.nativeToken());
        end
        function delete(obj)
            if obj.Token~=0
                radia_cuda_hmatrix('destroy_bh_material',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
    methods (Access=?radia.gpu.HMatrixOperator)
        function token=nativeToken(obj)
            assert(obj.Token~=0,'radia:gpu:BHMaterial','Material has been closed');
            token=obj.Token;
        end
    end
end
