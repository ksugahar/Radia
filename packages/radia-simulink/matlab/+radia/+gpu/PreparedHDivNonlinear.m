classdef PreparedHDivNonlinear < handle
    % Persistent shared-native HDiv BH/Newton workflow; NGSolve owns FE setup.
    properties (Access=private)
        State = struct()
    end
    methods
        function obj=PreparedHDivNonlinear(meshPath,order,appliedField,bhTable,options)
            arguments
                meshPath (1,:) char
                order (1,1) double {mustBeMember(order,[1 2])}
                appliedField (1,3) double {mustBeReal,mustBeFinite}
                bhTable (:,2) double {mustBeReal,mustBeFinite}
                options.Assembly (1,:) char {mustBeMember(options.Assembly,{'native-tet','ngsolve-host-normalized'})}='native-tet'
                options.HostSetup (1,1) struct=struct()
                options.NonlinearTolerance (1,1) double {mustBePositive,mustBeFinite}=5e-14
                options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
                options.LinearTolerance (1,1) double {mustBePositive,mustBeFinite}=5e-14
                options.LinearMaxIterations (1,1) double {mustBeInteger,mustBePositive}=10000
                options.MaterialBudgetBytes (1,1) double {mustBePositive,mustBeInteger}=512*2^20
                options.CacheBytes (1,1) double {mustBePositive,mustBeInteger}=32*2^20
                options.MaterialChunkElements (1,1) double {mustBeInteger,mustBeNonnegative}=0
                options.InexactNewton (1,1) logical=false
                options.FinalResidualOnly (1,1) logical=false
                options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'reference-diagonal','reference-chebyshev','reference-chebyshev-graph'})}='reference-diagonal'
                options.ChebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=32
                options.ChebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}=.01
            end
            assert(options.MaterialChunkElements==0 && ~options.InexactNewton && ~options.FinalResidualOnly, ...
                'radia:gpu:SharedNewtonOptions','Shared Newton does not implement chunked staging, inexact Newton or residual-only scheduling');
            assert(options.ChebyshevDegree<=64 && options.ChebyshevLowerFraction<1 && ...
                options.NonlinearTolerance<=1e-6 && options.LinearTolerance<=1e-6, ...
                'radia:gpu:SharedNewtonOptions','Invalid shared Newton accuracy or polynomial settings');
            normalized=strcmp(options.Assembly,'ngsolve-host-normalized');
            if normalized
                assert(isempty(meshPath),'radia:gpu:HostSetup','HostSetup owns checked mesh identity; meshPath must be empty');
                [setup,q,config]=radia.gpu.PreparedHDivNonlinear.readHostSetup(options.HostSetup,order,appliedField);
                meshHash=options.HostSetup.mesh_sha256;
            else
                assert(isempty(fieldnames(options.HostSetup)),'radia:gpu:HostSetup','HostSetup requires explicit ngsolve-host-normalized assembly');
                meshHash='';
            end
            % FE payload was checked above. Native leaf/key/charge-map checks
            % run after GPU selection but before resident allocation.
            whole=tic;device=gpuDevice;reserve=max(2*2^30,.2*device.TotalMemory);
            phase=tic;law=radia.gpu.makeHDivBHMaterial(bhTable);
            if ~normalized
                setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,order,appliedField);
                q=radia_cuda_projection('hdiv_material_quadrature',meshPath,order);
                charges=numel(setup.host);
            else
                charges=config.dimension;
            end
            nd=numel(setup.rhs);nq=numel(q.weights);
            assert(nd<=4096 && nq<=65536 && charges<=4608,'radia:gpu:SharedNewtonBudget','Shared FE staging bound exceeded');
            materialBound=32*nnz(q.evaluation)+256*nq+256*nd+32*2^20;
            combinedBound=materialBound+16*charges*charges+2*options.CacheBytes+128*2^20;
            assert(options.MaterialBudgetBytes<=512*2^20 && options.CacheBytes<=2*2^30 && ...
                materialBound<=options.MaterialBudgetBytes && combinedBound<=2*2^30 && ...
                device.AvailableMemory>reserve+combinedBound,'radia:gpu:SharedNewtonBudget','Shared resident array bound exceeded');
            if ~normalized
                if setup.curved_cells>0
                    cfg=radia.gpu.curvedHDivQuadrature(order);cfg.geometry_mode='curved';
                else
                    cfg=radia.gpu.affineHDivQuadrature(order);cfg.geometry_mode='affine';
                end
                cfg.schema_version=1;
                for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
                    cfg.(name{1})=setup.(name{1});
                end
                config=struct('assembly','hdiv_hacapk_device_dense','operator',cfg,'epsilon',1e-10, ...
                    'leaf_size',64,'eta',2.,'max_rank',64,'cache_bytes',options.CacheBytes);
            end
            feSeconds=toc(phase);phase=tic;
            op=radia.gpu.HMatrixOperator(config);material=[];
            try
                geometrySeconds=toc(phase);phase=tic;
                material=radia.gpu.HDivBHMaterial(q.evaluation,q.weights,law);
                wait(device);materialSeconds=toc(phase);
                obj.State=struct('op',op,'material',material,'mass',setup.mass,'rhs',setup.rhs, ...
                    'mean_operator',setup.mean_operator,'nd',nd,'nq',nq,'quadratureOrder',q.integration_order, ...
                    'barrierCap',law.cap,'options',rmfield(options,'HostSetup'),'assembly',options.Assembly,'meshHash',meshHash,'device',device,'combinedBound',combinedBound, ...
                    'materialBound',materialBound,'prepare_seconds',toc(whole), ...
                    'prepare_phase_seconds',struct('fe_and_law',feSeconds,'geometry_create',geometrySeconds,'material_upload',materialSeconds));
            catch ex
                if ~isempty(material),delete(material);end
                delete(op);rethrow(ex);
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
            native=s.op.solveNonlinearMaterial(s.material,s.mass,scale*s.rhs,initialCoefficients,config);
            value=s.material.evaluate(native.x);wait(s.device);
            result=struct('coefficients',native.x,'mean_magnetization',s.mean_operator*native.x, ...
                'nonlinear_relative_residual',native.relative_residual,'nonlinear_iterations',native.newton_iterations, ...
                'history',native.history,'final_energy',native.energy,'maximum_material_magnetization',value.maximum_magnetization, ...
                'material_barrier_cap',s.barrierCap,'operator_stats',s.op.stats(),'material_stats',s.material.stats(), ...
                'preconditioner_graph_launches',native.preconditioner_graph_launches,'preconditioner',o.Preconditioner, ...
                'whole_seconds',toc(whole),'prepare_seconds',s.prepare_seconds,'prepare_phase_seconds',s.prepare_phase_seconds, ...
                'prepared_reuse',true,'field_scale',scale,'quadrature_order',s.quadratureOrder,'material_points',s.nq, ...
                'combined_array_bound_bytes',s.combinedBound,'material_array_bound_bytes',s.materialBound, ...
                'nonlinear',true,'backend','gpu','implementation','shared-native','assembly_route',s.assembly,'mesh_sha256',s.meshHash, ...
                'gram_integration_profile','Host-normalized follows declared Radia charge builder; native-tet uses all-high-order; independent CPU reference must match the selected assembly', ...
                'scope','Persistent shared BH and frozen-tangent Newton; solve excludes preparation/destruction', ...
                'transfers','Fixed FE/weights/law upload once; shared native Newton retains device state and workspaces', ...
                'memory_limitations','Conservative array bound and native owner counters; not host/library/global peak');
        end
        function info=stats(obj)
            assert(~isempty(fieldnames(obj.State)),'Prepared object has been released');
            s=obj.State;info=struct('prepare_seconds',s.prepare_seconds,'prepare_phase_seconds',s.prepare_phase_seconds, ...
                'combined_array_bound_bytes',s.combinedBound,'operator',s.op.stats(),'material',s.material.stats(),'implementation','shared-native');
        end
        function delete(obj)
            if isfield(obj.State,'op')
                s=obj.State;obj.State=struct();
                try,delete(s.material);finally,delete(s.op);end
            end
        end
    end
    methods (Static,Access=private)
        function [setup,q,config]=readHostSetup(h,order,field)
            required={'schema';'backend';'order';'mesh_sha256';'label_contract_sha256';'operator';'mass';'rhs_operator';'mean_operator';'evaluation';'weights';'quadrature_order'};
            assert(isequal(sort(fieldnames(h)),sort(required)) && strcmp(h.schema,'radia.hdiv-host-setup.v1') && ...
                strcmp(h.backend,'gpu') && isequal(h.order,order),'radia:gpu:HostSetup','Invalid explicit host setup schema/order/backend');
            assert(ischar(h.mesh_sha256) && ~isempty(regexp(h.mesh_sha256,'^[0-9a-f]{64}$','once')) && ...
                ischar(h.label_contract_sha256) && ~isempty(regexp(h.label_contract_sha256,'^[0-9a-f]{64}$','once')), ...
                'radia:gpu:HostSetup','Missing checked mesh/label identity');
            nd=size(h.mass,1);nq=numel(h.weights);
            assert(nd>0 && nd<=4096 && size(h.mass,2)==nd && issparse(h.mass) && isa(h.mass,'double') && ...
                isreal(h.mass) && all(isfinite(nonzeros(h.mass))) && nnz(h.mass)<=4000000, ...
                'radia:gpu:HostSetup','Invalid sparse NGSolve mass');
            assert(isa(h.rhs_operator,'double') && isreal(h.rhs_operator) && isequal(size(h.rhs_operator),[nd 3]) && ...
                all(isfinite(h.rhs_operator(:))) && isa(h.mean_operator,'double') && isreal(h.mean_operator) && ...
                isequal(size(h.mean_operator),[3 nd]) && all(isfinite(h.mean_operator(:))), ...
                'radia:gpu:HostSetup','Invalid NGSolve excitation/mean numbering');
            assert(nq>0 && nq<=65536 && isa(h.weights,'double') && iscolumn(h.weights) && isreal(h.weights) && ...
                all(isfinite(h.weights) & h.weights>0) && issparse(h.evaluation) && isa(h.evaluation,'double') && ...
                isreal(h.evaluation) && isequal(size(h.evaluation),[3*nq nd]) && nnz(h.evaluation)<=4000000 && ...
                all(isfinite(nonzeros(h.evaluation))) && isequal(h.quadrature_order,6.), ...
                'radia:gpu:HostSetup','Invalid NGSolve mapped quadrature');
            config=h.operator;
            assert(isstruct(config) && isscalar(config) && isequal(config.schema_version,2.) && ...
                strcmp(config.assembly,'normalized_leaves') && strcmp(config.operator_kind,'demag') && ...
                config.dimension>0 && config.dimension<=4608 && issparse(config.charge_map) && ...
                isequal(size(config.charge_map),[config.dimension nd]), ...
                'radia:gpu:HostSetup','Invalid normalized charge operator');
            setup=struct('mass',h.mass,'rhs',h.rhs_operator*field(:),'mean_operator',h.mean_operator);
            q=struct('evaluation',h.evaluation,'weights',h.weights,'integration_order',h.quadrature_order);
        end
    end

end
