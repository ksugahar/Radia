classdef HDivDeviceMass < handle
    %HDIVDEVICEMASS Device-built reference mass E'*W*E of a device setup.
    % The NGSolve material quadrature tables build M on the device; the
    % linear material solve uses it in place, with no host mass upload.
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj=HDivDeviceMass(setup)
            arguments
                setup (1,1) struct
            end
            obj.Token=radia_hdiv_device('create_mass',setup.evaluation,setup.dofs);
        end
        function result=apply(obj,x)
            arguments
                obj (1,1) radia.gpu.HDivDeviceMass
                x (:,1) double {mustBeReal,mustBeFinite}
            end
            result=radia_hdiv_device('mass_apply',obj.nativeToken(),x);
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
            assert(obj.Token~=0,'radia:gpu:HDivDevice','Mass has been closed');
            token=obj.Token;
        end
    end
end
