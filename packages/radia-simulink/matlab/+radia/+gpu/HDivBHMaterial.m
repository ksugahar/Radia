classdef HDivBHMaterial < handle
    %HDIVBHMATERIAL Persistent device inverse-BH material of a device setup.
    % The evaluation map and reference mass are built on the device from the
    % NGSolve quadrature tables; only states, directions and results cross.
    % This owner does not assemble FE geometry or select a CPU fallback.
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj=HDivBHMaterial(setup,law)
            arguments
                setup (1,1) struct
                law (1,1) struct
            end
            obj.Token=radia_hdiv_device('create_material',setup.evaluation,setup.dofs,law);
        end
        function result=evaluate(obj,state,direction)
            % Energy, gradient and maximum |M|; with a direction also the tangent action.
            arguments
                obj (1,1) radia.gpu.HDivBHMaterial
                state (:,1) double {mustBeReal,mustBeFinite}
                direction (:,1) double {mustBeReal,mustBeFinite}=[]
            end
            result=radia_hdiv_device('material_evaluate',obj.nativeToken(),state,direction);
        end
        function info=stats(obj)
            info=radia_hdiv_device('material_stats',obj.nativeToken());
        end
        function delete(obj)
            if obj.Token~=0
                radia_hdiv_device('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
    methods (Access=?radia.gpu.HDivDeviceOperator)
        function token=nativeToken(obj)
            assert(obj.Token~=0,'radia:gpu:BHMaterial','Material has been closed');
            token=obj.Token;
        end
    end
end
