classdef PreparedHDivNonlinear < handle
    % Persistent device HDiv BH/Newton workflow over a radia.hdiv-device-setup.v1.
    % NGSolve owns basis, orientation and quadrature (staged by Radia Python);
    % the Gram, charge map, evaluation map and reference mass are built on the
    % device. Only vectors and scalar results cross per solve.
    properties (Access=private)
        State = struct()
    end
    methods
        function obj=PreparedHDivNonlinear(setup,appliedField,bhTable,options)
            arguments
                setup (1,1) struct
                appliedField (1,3) double {mustBeReal,mustBeFinite}
                bhTable (:,2) double {mustBeReal,mustBeFinite}
                options.NonlinearTolerance (1,1) double {mustBePositive,mustBeFinite}=5e-14
                options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
                options.LinearTolerance (1,1) double {mustBePositive,mustBeFinite}=5e-14
                options.LinearMaxIterations (1,1) double {mustBeInteger,mustBePositive}=10000
                options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'reference-diagonal','reference-chebyshev','reference-chebyshev-graph'})}='reference-diagonal'
                options.ChebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=32
                options.ChebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}=.01
            end
            assert(options.ChebyshevDegree<=64 && options.ChebyshevLowerFraction<1 && ...
                options.NonlinearTolerance<=1e-6 && options.LinearTolerance<=1e-6, ...
                'radia:gpu:SharedNewtonOptions','Invalid shared Newton accuracy or polynomial settings');
            setup=radia.gpu.readHDivDeviceSetup(setup);
            whole=tic;device=gpuDevice;
            phase=tic;law=radia.gpu.makeHDivBHMaterial(bhTable);lawSeconds=toc(phase);
            % Device budgets are checked natively before any resident allocation.
            phase=tic;material=radia.gpu.HDivBHMaterial(setup,law);op=[];
            try
                wait(device);materialSeconds=toc(phase);phase=tic;
                op=radia.gpu.HDivDeviceOperator(setup);wait(device);operatorSeconds=toc(phase);
                obj.State=struct('op',op,'material',material,'rhs',setup.rhs_operator*appliedField(:), ...
                    'mean_operator',setup.mean_operator,'nd',setup.dofs,'order',setup.order,'barrierCap',law.cap, ...
                    'options',options,'meshHash',setup.mesh_sha256,'device',device,'assembly',op.assembly(), ...
                    'prepare_seconds',toc(whole),'prepare_phase_seconds', ...
                    struct('law',lawSeconds,'material_create',materialSeconds,'operator_create',operatorSeconds));
            catch ex
                if ~isempty(op),delete(op);end
                delete(material);rethrow(ex);
            end
        end
        function result=solve(obj,scale,initialCoefficients)
            arguments
                obj (1,1) radia.gpu.PreparedHDivNonlinear
                scale (1,1) double {mustBeReal,mustBeFinite}=1
                initialCoefficients (:,1) double {mustBeReal,mustBeFinite}=[]
            end
            assert(~isempty(fieldnames(obj.State)),'Prepared object has been released');
            whole=tic;s=obj.State;o=s.options;
            if isempty(initialCoefficients),initialCoefficients=zeros(s.nd,1);end
            assert(numel(initialCoefficients)==s.nd,'Invalid explicit initial coefficients');
            config=struct('nonlinear_tolerance',o.NonlinearTolerance,'linear_tolerance',o.LinearTolerance, ...
                'maximum_iterations',o.MaxIterations,'maximum_linear_iterations',o.LinearMaxIterations, ...
                'preconditioner',o.Preconditioner,'degree',o.ChebyshevDegree,'lower_fraction',o.ChebyshevLowerFraction);
            native=s.op.solveNonlinear(s.material,scale*s.rhs,initialCoefficients,config);
            value=s.material.evaluate(native.x);wait(s.device);
            result=struct('coefficients',native.x,'mean_magnetization',s.mean_operator*native.x, ...
                'nonlinear_relative_residual',native.relative_residual,'nonlinear_iterations',native.newton_iterations, ...
                'history',native.history,'final_energy',native.energy,'maximum_material_magnetization',value.maximum_magnetization, ...
                'material_barrier_cap',s.barrierCap,'operator_stats',s.op.stats(),'material_stats',s.material.stats(), ...
                'preconditioner_graph_launches',native.preconditioner_graph_launches,'preconditioner',o.Preconditioner, ...
                'whole_seconds',toc(whole),'prepare_seconds',s.prepare_seconds,'prepare_phase_seconds',s.prepare_phase_seconds, ...
                'prepared_reuse',true,'field_scale',scale,'order',s.order,'dofs',s.nd, ...
                'nonlinear',true,'backend','gpu','assembly_route',s.assembly.route,'mesh_sha256',s.meshHash, ...
                'scope','Persistent device BH and frozen-tangent Newton; solve excludes preparation/destruction', ...
                'transfers','Primitive NGSolve tables and the law upload once; per solve only vectors cross');
        end
        function info=stats(obj)
            assert(~isempty(fieldnames(obj.State)),'Prepared object has been released');
            s=obj.State;info=struct('prepare_seconds',s.prepare_seconds,'prepare_phase_seconds',s.prepare_phase_seconds, ...
                'assembly',s.assembly,'operator',s.op.stats(),'material',s.material.stats());
        end
        function delete(obj)
            if isfield(obj.State,'op')
                s=obj.State;obj.State=struct();
                try,delete(s.material);finally,delete(s.op);end
            end
        end
    end
end
