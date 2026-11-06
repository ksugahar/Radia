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
        function result=linearize(obj,state)
            % Freeze the tangent at state; result.linearization_epoch names it.
            arguments
                obj (1,1) radia.gpu.HDivBHMaterial
                state (:,1) double {mustBeReal,mustBeFinite}
            end
            result=radia_hdiv_device('material_linearize',obj.nativeToken(),state);
        end
        function action=applyTangent(obj,direction,epoch)
            % Tangent action of a current linearization; a stale epoch is refused.
            arguments
                obj (1,1) radia.gpu.HDivBHMaterial
                direction (:,1) double {mustBeReal,mustBeFinite}
                epoch (1,1) double {mustBeInteger,mustBePositive}
            end
            action=radia_hdiv_device('material_apply_tangent',obj.nativeToken(),direction,epoch);
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
    methods (Static)
        function result=inverseSamples(law,amplitudes)
            % Device inverse-BH secant, tangent and energy at the given amplitudes.
            arguments
                law (1,1) struct
                amplitudes (:,1) double {mustBeReal,mustBeFinite}
            end
            result=radia_hdiv_device('inverse_bh_samples',law,amplitudes);
        end
    end
    methods (Access=?radia.gpu.HDivDeviceOperator)
        function token=nativeToken(obj)
            assert(obj.Token~=0,'radia:gpu:BHMaterial','Material has been closed');
            token=obj.Token;
        end
    end
end
