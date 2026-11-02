classdef HDivDeviceOperator < handle
    %HDIVDEVICEOPERATOR Device-built charge Gram and charge map B'*G*B.
    % Built on the device from a radia.hdiv-device-setup.v1 payload; only
    % vectors and scalar results cross. Solves use device-owned masses and
    % materials of the same setup; no matrix is uploaded from MATLAB.
    properties (Access=private)
        Token = uint64(0)
    end
    properties (SetAccess=private)
        Dofs (1,1) double = 0
    end
    methods
        function obj=HDivDeviceOperator(setup)
            arguments
                setup (1,1) struct
            end
            obj.Token=radia_hdiv_device('create_operator',setup.charges);
            obj.Dofs=setup.dofs;
        end
        function y=apply(obj,x)
            arguments
                obj (1,1) radia.gpu.HDivDeviceOperator
                x (:,1) double {mustBeReal,mustBeFinite}
            end
            y=radia_hdiv_device('apply',obj.nativeToken(),x);
        end
        function info=assembly(obj)
            info=radia_hdiv_device('assembly',obj.nativeToken());
        end
        function info=stats(obj)
            info=radia_hdiv_device('stats',obj.nativeToken());
        end
        function result=solveLinear(obj,mass,rhs,invChi,tolerance,maxIterations,preconditioner,chebyshevDegree,chebyshevLowerFraction)
            arguments
                obj (1,1) radia.gpu.HDivDeviceOperator
                mass (1,1) radia.gpu.HDivDeviceMass
                rhs (:,1) double {mustBeReal,mustBeFinite}
                invChi (1,1) double {mustBePositive,mustBeFinite}
                tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
                maxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
                preconditioner (1,:) char {mustBeMember(preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial','mass-polynomial-graph','mass-chebyshev','mass-chebyshev-graph'})}='mass-polynomial'
                chebyshevDegree double {mustBeScalarOrEmpty,mustBeInteger,mustBePositive}=[]
                chebyshevLowerFraction double {mustBeScalarOrEmpty,mustBePositive,mustBeFinite}=[]
            end
            % Chebyshev settings are accepted only with a Chebyshev mode (native check).
            extra={};
            if ~isempty(chebyshevDegree) || ~isempty(chebyshevLowerFraction)
                assert(~isempty(chebyshevDegree) && ~isempty(chebyshevLowerFraction), ...
                    'radia:gpu:HDivLinear','Chebyshev degree and lower fraction are given together');
                extra={chebyshevDegree,chebyshevLowerFraction};
            end
            result=radia_hdiv_device('solve_linear',obj.nativeToken(),mass.nativeToken(),rhs,invChi,tolerance,maxIterations,preconditioner,extra{:});
        end
        function result=solveNonlinear(obj,material,rhs,initial,options)
            % Shared frozen-tangent Newton and line search over the material's device mass.
            arguments
                obj (1,1) radia.gpu.HDivDeviceOperator
                material (1,1) radia.gpu.HDivBHMaterial
                rhs (:,1) double {mustBeReal,mustBeFinite}
                initial (:,1) double {mustBeReal,mustBeFinite}
                options (1,1) struct
            end
            result=radia_hdiv_device('solve_nonlinear',obj.nativeToken(),material.nativeToken(),rhs,initial,options);
        end
        function delete(obj)
            if obj.Token~=0
                radia_hdiv_device('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
    methods (Access=private)
        function token=nativeToken(obj)
            assert(obj.Token~=0,'radia:gpu:HDivDevice','Operator has been closed');
            token=obj.Token;
        end
    end
end
