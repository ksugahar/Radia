classdef PreparedHDivNonlinear < handle
    % Explicit fixed-mesh, fixed-BH owner; no hidden cache or accepted-state mutation.
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
    options.NonlinearTolerance (1,1) double {mustBePositive}=1e-6
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
    options.LinearTolerance (1,1) double {mustBePositive}=1e-10
    options.LinearMaxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
    options.MaterialBudgetBytes (1,1) double {mustBePositive,mustBeInteger}=512*2^20
    options.CacheBytes (1,1) double {mustBePositive,mustBeInteger}=32*2^20
    options.MaterialChunkElements (1,1) double {mustBeInteger,mustBeNonnegative}=0
    options.InexactNewton (1,1) logical=false
    options.FinalResidualOnly (1,1) logical=false
    options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial','mass-polynomial-graph','mass-chebyshev','mass-chebyshev-graph'})}='mass-cholesky'
    options.ChebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=8
    options.ChebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}=.1
end
            obj.State=prepareState(meshPath,order,appliedField,bhTable,options);
        end
        function result=solve(obj,scale,initialCoefficients)
            arguments
                obj (1,1) radia.gpu.PreparedHDivNonlinear
                scale (1,1) double {mustBeReal,mustBeFinite}=1
                initialCoefficients (:,1) double {mustBeReal,mustBeFinite}=[]
            end
            assert(~isempty(fieldnames(obj.State)),'Prepared object has been released');
            result=solveState(obj.State,scale,initialCoefficients);
        end
        function info=stats(obj)
            assert(~isempty(fieldnames(obj.State)),'Prepared object has been released');
            info=struct('prepare_seconds',obj.State.prepare_seconds,'prepare_phase_seconds',obj.State.prepare_phase_seconds,'combined_array_bound_bytes',obj.State.combinedBound,'operator',obj.State.op.stats());
        end
        function delete(obj)
            if isfield(obj.State,'op')
                native=obj.State.op;obj.State=struct();delete(native);
            end
        end
    end
end

function state=prepareState(meshPath,order,appliedField,bhTable,options)
assert(options.ChebyshevDegree<=64 && options.ChebyshevLowerFraction<1,'Invalid fixed Chebyshev settings');
assert(ismember(options.Preconditioner,{'mass-chebyshev','mass-chebyshev-graph'}) || (options.ChebyshevDegree==8 && options.ChebyshevLowerFraction==.1),'Chebyshev settings require mass-chebyshev');
assert(options.NonlinearTolerance<=1e-6 && options.LinearTolerance<=1e-6);
whole=tic;device=gpuDevice;reserve=max(2*2^30,.2*device.TotalMemory);
assert(device.AvailableMemory>reserve+256*2^20,'Nonlinear GPU reserve');deviceInitSeconds=toc(whole);phase=tic;
law=radia.gpu.makeHDivBHMaterial(bhTable);barrierCap=law.cap;
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,order,appliedField);
nd=numel(setup.rhs);chunks={};E={};weights=[];hostChunks={};
% Native cache capacity is 2 GiB; the combined phase bound below still
% rejects requests whose cache plus resident operator exceed 2 GiB.
assert(options.MaterialBudgetBytes<=2^30 && options.CacheBytes<=2*2^30 && options.MaterialChunkElements<=128,'Explicit material/cache/chunk budgets');
if options.MaterialChunkElements==0
    q=radia_cuda_projection('hdiv_material_quadrature',meshPath,order);nq=numel(q.weights);quadratureOrder=q.integration_order;
    assert(nd<=4096 && nq<=65536,'Unchunked nonlinear staging limit');
    payload=24*nnz(q.evaluation)+8*(3*nq+9*nd*nd)+32*2^20;
    combinedBound=payload+2*options.CacheBytes+32*nd*nd+64*2^20;
    chunkCount=1;shapeEntries=nnz(q.evaluation);tangentBound=nd*nd;
else
    assert(nd<=16384,'Chunked nonlinear DOF limit');ne=size(setup.cell_nodes,1);nq=0;shapeEntries=0;tangentBound=0;maxChunkPoints=0;seen=[];
    materialContext=radia_cuda_projection('create_hdiv_material_context',meshPath,order);
    materialOwner=onCleanup(@()radia_cuda_projection('destroy_hdiv_material_context',materialContext));
    for first=0:options.MaterialChunkElements:ne-1
        part=radia_cuda_projection('hdiv_material_context_chunk',materialContext,first,min(options.MaterialChunkElements,ne-first));
        if isempty(seen),assert(part.total_points<=2000000);seen=false(part.total_points,1);quadratureOrder=part.integration_order;end
        ids=part.global_point_indices(:)+1;assert(part.dofs==nd && part.total_points==numel(seen) && part.integration_order==quadratureOrder && all(ids>=1 & ids<=numel(seen)) && numel(unique(ids))==numel(ids) && ~any(seen(ids)));seen(ids)=true;
        nq=nq+numel(ids);maxChunkPoints=max(maxChunkPoints,numel(ids));shapeEntries=shapeEntries+nnz(part.evaluation);tangentBound=tangentBound+part.tangent_nonzero_bound;
        % Stop before retaining more chunks when even this partial lower bound
        % exceeds the requested final array allowance. Native staging still
        % bounds the current chunk; this does not claim a host/global peak.
        partialPayload=24*shapeEntries+288*min(tangentBound,nd*nd)+384*maxChunkPoints+64*2^20;
        assert(partialPayload<=options.MaterialBudgetBytes,'radia:gpu:MaterialChunkBudget','Material chunk array budget exceeded during staging');
        hostChunks{end+1}=part;
    end
    assert(all(seen));clear materialOwner materialContext seen part
    tangentBound=min(tangentBound,nd*nd);chunkCount=numel(hostChunks);
    % Every E_i' D E_j couples only DOFs sharing a staged element.
    % Twelve 24-byte sparse-array allowances plus bounded point temporaries;
    % CUDA library workspace and global peak are not covered by this array bound.
    payload=24*shapeEntries+288*tangentBound+384*maxChunkPoints+64*2^20;
    charges=numel(setup.host);assert(charges<=4608,'Resident geometry source budget');
    residentBound=16*charges*charges+64*2^20;
    if strcmp(options.Preconditioner,'mass-cholesky'),solverBound=32*nd*nd;else,solverBound=128*tangentBound+512*nd+64*2^20;end
    % Deferred oracle buffers are released inside MEX before handle publication.
    % Material GPU maps are uploaded only after that geometry setup completes.
    geometryBound=options.CacheBytes+residentBound;
    materialSolveBound=payload+residentBound+solverBound;
    combinedBound=max(geometryBound,materialSolveBound);
end
assert(payload<=options.MaterialBudgetBytes && combinedBound<=2*2^30 && device.AvailableMemory>reserve+combinedBound,'Nonlinear combined array budget');
if setup.curved_cells>0,cfg=radia.gpu.curvedHDivQuadrature(order);cfg.geometry_mode='curved';else,cfg=radia.gpu.affineHDivQuadrature(order);cfg.geometry_mode='affine';end
cfg.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'},cfg.(name{1})=setup.(name{1});end
config=struct('assembly','hdiv_hacapk_device_dense','operator',cfg,'epsilon',1e-10,'leaf_size',64,'eta',2.,'max_rank',64,'cache_bytes',options.CacheBytes);
feSeconds=toc(phase);phase=tic;op=radia.gpu.HMatrixOperator(config);geometrySeconds=toc(phase);phase=tic;
try
if options.MaterialChunkElements==0
    nq=numel(q.weights);E=cell(1,3);for c=1:3,E{c}=gpuArray(q.evaluation((c-1)*nq+(1:nq),:));end
    weights=gpuArray(q.weights);
else
    for k=1:numel(hostChunks)
        part=hostChunks{k};count=numel(part.weights);item=struct('E',{{}},'weights',gpuArray(part.weights));
        for c=1:3,item.E{c}=gpuArray(part.evaluation((c-1)*count+(1:count),:));end
        chunks{end+1}=item;
    end
    clear item part hostChunks
end
law.device=struct('knots',gpuArray(law.knots),'coefficients',gpuArray(law.coefficients),'integral_base',gpuArray(law.integral_base));
lawDeviceBytes=8*(numel(law.knots)+numel(law.coefficients)+numel(law.integral_base));
uploadSeconds=toc(phase);wait(device);
catch ex
    delete(op);rethrow(ex);
end
state=struct();
state.device=device;
state.law=law;
state.barrierCap=barrierCap;
state.setup=setup;
state.nd=nd;
state.chunks=chunks;
state.E=E;
state.weights=weights;
state.op=op;
state.nq=nq;
state.quadratureOrder=quadratureOrder;
state.chunkCount=chunkCount;
state.shapeEntries=shapeEntries;
state.tangentBound=tangentBound;
state.payload=payload;
state.combinedBound=combinedBound;
state.lawDeviceBytes=lawDeviceBytes;
state.prepare_seconds=toc(whole);state.options=options;state.prepare_phase_seconds=struct('device_initialization',deviceInitSeconds,'fe_and_law',feSeconds,'geometry_create',geometrySeconds,'material_upload_enqueue',uploadSeconds);

end

function result=solveState(state,scale,initialCoefficients)
whole=tic;
device=state.device;
law=state.law;
barrierCap=state.barrierCap;
setup=state.setup;
nd=state.nd;
chunks=state.chunks;
E=state.E;
weights=state.weights;
op=state.op;
nq=state.nq;
quadratureOrder=state.quadratureOrder;
chunkCount=state.chunkCount;
shapeEntries=state.shapeEntries;
tangentBound=state.tangentBound;
payload=state.payload;
combinedBound=state.combinedBound;
lawDeviceBytes=state.lawDeviceBytes;
options=state.options;reserve=max(2*2^30,.2*device.TotalMemory);assert(device.AvailableMemory>reserve+combinedBound,'Prepared nonlinear transient reserve');
phase=tic;rhs=scale*setup.rhs;sourceNorm=norm(rhs);m=zeros(nd,1);maxSolverOwned=0;graphLaunches=0;warmIterations=0;warmTrueResidual=0;warmTolerance=options.LinearTolerance;if options.InexactNewton,warmTolerance=max(warmTolerance,min(1e-6,sqrt(options.NonlinearTolerance)));end
if ~isempty(initialCoefficients)
 assert(numel(initialCoefficients)==nd && all(isfinite(initialCoefficients)),'Invalid explicit initial coefficients');m=initialCoefficients(:);
else
if sourceNorm>0
    if isempty(chunks)
        mass0=E{1}'*spdiags(weights,0,nq,nq)*E{1}+E{2}'*spdiags(weights,0,nq,nq)*E{2}+E{3}'*spdiags(weights,0,nq,nq)*E{3};
        mass0=triu(mass0)+triu(mass0,1)';
        initial=op.solveMaterial(gather(mass0),rhs,law.initial_reluctivity,warmTolerance,options.LinearMaxIterations,options.Preconditioner,options.ChebyshevDegree,options.ChebyshevLowerFraction);clear mass0
    else
        initialMaterial=radia.gpu.evaluateHDivBHChunks(chunks,law,m,true);
        initial=op.solveMaterial(initialMaterial.tangent,rhs,1,warmTolerance,options.LinearMaxIterations,options.Preconditioner,options.ChebyshevDegree,options.ChebyshevLowerFraction);clear initialMaterial
    end
    m=initial.x;maxSolverOwned=initial.solver_owned_device_bytes;warmIterations=initial.iterations;warmTrueResidual=initial.true_relative_residual;
    if isfield(initial,'preconditioner_graph_launches'),graphLaunches=graphLaunches+initial.preconditioner_graph_launches;end
end
end
warmSeconds=toc(phase);materialSeconds=0;linearSeconds=0;lineSeconds=0;
history=struct('relative_residual',{},'energy',{},'alpha',{},'backtracks',{},'inner_iterations',{},'inner_true_residual',{},'inner_tolerance',{},'full_trial_relative_residual',{},'full_trial_energy_gap_ulps',{},'accepted_trial_relative_residual',{},'acceptance_reason',{});
converged=false;certifiedTrial=false;outerTangents=0;residualOnlyChecks=0;
for iteration=0:options.MaxIterations
    wantTangent=~(options.FinalResidualOnly && certifiedTrial);
    phase=tic;mat=materialValue(chunks,E,weights,law,m,wantTangent);outerTangents=outerTangents+wantTangent;residualOnlyChecks=residualOnlyChecks+~wantTangent;
    op.upload(m);op.apply();nm=op.download();R=mat.gradient+nm-rhs;
    relative=norm(R)/max(sourceNorm,realmin);energy=mat.energy+.5*m'*nm-rhs'*m;
    if ~wantTangent && relative>options.NonlinearTolerance
        % Never reuse the trial certificate as convergence proof. If this
        % independent residual check fails, construct the tangent to continue.
        mat=materialValue(chunks,E,weights,law,m,true);outerTangents=outerTangents+1;
        R=mat.gradient+nm-rhs;relative=norm(R)/max(sourceNorm,realmin);energy=mat.energy+.5*m'*nm-rhs'*m;
    end
    materialSeconds=materialSeconds+toc(phase);
    assert(isfinite(relative) && isfinite(energy));
    if relative<=options.NonlinearTolerance,converged=true;break,end
    if iteration>=options.MaxIterations
        diagnostic=struct('schema','radia.nonlinear-convergence.v1','iteration',iteration,'current_relative_residual',relative,'target',options.NonlinearTolerance,'energy',energy,'state_norm',norm(m),'history',history);
        error('radia:gpu:NewtonConvergence','Nonlinear assembled residual failed convergence; diagnostics=%s',jsonencode(diagnostic));
    end
    innerTolerance=options.LinearTolerance;if options.InexactNewton,innerTolerance=max(innerTolerance,min(1e-6,.1*relative));end
    phase=tic;linear=op.solveMaterial(mat.tangent,-R,1,innerTolerance,options.LinearMaxIterations,options.Preconditioner,options.ChebyshevDegree,options.ChebyshevLowerFraction);linearSeconds=linearSeconds+toc(phase);
    if isfield(linear,'preconditioner_graph_launches'),graphLaunches=graphLaunches+linear.preconditioner_graph_launches;end
    maxSolverOwned=max(maxSolverOwned,linear.solver_owned_device_bytes);delta=linear.x;decrement=-R'*delta;
    assert(isfinite(decrement) && decrement>0,'Newton direction is not descending');
    phase=tic;alpha=1;backtracks=0;accepted=false;trialTrace=zeros(0,5);acceptanceReason='armijo';
    while alpha>1e-10
        candidate=m+alpha*delta;
        trial=materialValue(chunks,E,weights,law,candidate,false);op.upload(candidate);op.apply();trialN=op.download();trialEnergy=trial.energy+.5*candidate'*trialN-rhs'*candidate;
        trialRelative=norm(trial.gradient+trialN-rhs)/max(sourceNorm,realmin);
        armijoBound=energy-1e-4*alpha*decrement;energyUlp=eps(max(abs(energy),1));
        trialTrace(end+1,:)=[alpha,trialRelative,trialEnergy,armijoBound,(trialEnergy-armijoBound)/energyUlp];
        % A root already meeting the requested original-equation criterion
        % needs no globalization. Absolute coenergy subtraction can reject
        % such a root by a few ULPs. The next outer iteration recertifies R.
        if isfinite(trialEnergy) && isfinite(trialRelative) && trialRelative<=options.NonlinearTolerance
            accepted=true;acceptanceReason='original-residual-certified';break
        end
        if isfinite(trialEnergy) && trialEnergy<=energy-1e-4*alpha*decrement,accepted=true;break,end
        alpha=alpha*.5;backtracks=backtracks+1;
    end
    lineSeconds=lineSeconds+toc(phase);
    if ~accepted
        diagnostic=struct('schema','radia.nonlinear-line-search.v1','iteration',iteration,'current_relative_residual',relative, ...
            'target',options.NonlinearTolerance,'energy',energy,'predicted_decrease',decrement,'energy_ulp',energyUlp, ...
            'step_norm',norm(delta),'state_norm',norm(m),'field_scale',scale,'preconditioner',options.Preconditioner,'chebyshev_degree',options.ChebyshevDegree,'chebyshev_lower_fraction',options.ChebyshevLowerFraction, ...
            'inner_true_residual',linear.true_relative_residual,'material_max_M',mat.maximum_magnetization, ...
            'trial_columns',{{'alpha','original_relative_residual','energy','armijo_bound','energy_gap_ulps'}},'trials',trialTrace);
        error('radia:gpu:NewtonLineSearch','Energy line search failed; accepted state unchanged; diagnostics=%s',jsonencode(diagnostic));
    end
    history(end+1)=struct('relative_residual',relative,'energy',energy,'alpha',alpha,'backtracks',backtracks,'inner_iterations',linear.iterations,'inner_true_residual',linear.true_relative_residual,'inner_tolerance',innerTolerance,'full_trial_relative_residual',trialTrace(1,2),'full_trial_energy_gap_ulps',trialTrace(1,5),'accepted_trial_relative_residual',trialTrace(end,2),'acceptance_reason',acceptanceReason);
    m=candidate;certifiedTrial=strcmp(acceptanceReason,'original-residual-certified');
end
phase=tic;assert(converged);stats=op.stats();mean=setup.mean_operator*m;
% Drop per-solve aliases; the prepared object retains immutable GPU owners.
clear chunks E weights law setup op
wait(device);cleanupSeconds=toc(phase);
hostPhases=struct('warmstart',warmSeconds,'residual_tangent_and_demag',materialSeconds,'newton_linear_solve',linearSeconds,'line_search',lineSeconds,'cleanup',cleanupSeconds);
result=struct('coefficients',m,'mean_magnetization',mean,'nonlinear_relative_residual',relative, ...
    'maximum_material_magnetization',mat.maximum_magnetization,'material_barrier_cap',barrierCap, ...
    'nonlinear_iterations',numel(history),'history',history,'final_energy',energy,'operator_stats',stats, ...
    'maximum_solver_owned_device_bytes',maxSolverOwned,'whole_seconds',toc(whole),'quadrature_order',quadratureOrder, ...
    'warmstart_iterations',warmIterations,'warmstart_true_residual',warmTrueResidual,'warmstart_tolerance',warmTolerance, ...
    'preconditioner_graph_launches',graphLaunches, ...
    'final_residual_only',options.FinalResidualOnly,'outer_tangent_evaluations',outerTangents,'residual_only_outer_evaluations',residualOnlyChecks, ...
    'preconditioner',options.Preconditioner,'chebyshev_degree',options.ChebyshevDegree,'chebyshev_lower_fraction',options.ChebyshevLowerFraction, ...
    'prepared_reuse',true,'inexact_newton',options.InexactNewton,'field_scale',scale,'prepare_seconds',state.prepare_seconds, ...
    'material_chunk_count',chunkCount,'material_points',nq,'material_shape_entries',shapeEntries,'material_tangent_nonzero_bound',tangentBound, ...
    'host_phase_seconds',hostPhases,'phase_timing_scope','Host elapsed API regions; asynchronous upload work can carry into warmstart; not CUDA kernel timings', ...
    'energy_units','A^2*m; common mu0 omitted, multiply integrated coenergy by mu0 for joules', ...
    'material_array_bound_bytes',payload,'combined_array_bound_bytes',combinedBound,'law_owned_device_bytes',lawDeviceBytes,'nonlinear',true,'live_count',radia_cuda_hmatrix('live_count'), ...
    'gram_integration_profile','all-high-order; independent CPU reference must explicitly select ho_far_factor=inf', ...
    'trial_acceptance','Requested original nonlinear residual checked before energy globalization; nonconverged trials retain Armijo; final outer residual recertified', ...
    'scope','Prepared fixed-geometry and fixed-BH context; solve time excludes preparation and destruction', ...
    'transfers','Fixed FE/weights/shape maps and geometry retained by explicit owner; chunk material sums stay on GPU until each complete gradient/tangent is gathered; inverse-law coefficient arrays uploaded once per solve and reused; one coefficient upload per material evaluation; coefficient vectors transferred for resident demag and inner solves; no Python GPU adapter', ...
    'memory_limitations','Native owned counters and conservative MATLAB array bound; no CUDA/library/host/global peak');
end


function out=materialValue(chunks,E,weights,law,m,wantTangent)
if isempty(chunks),out=radia.gpu.evaluateHDivBHMaterial(E,weights,law,m,wantTangent);else,out=radia.gpu.evaluateHDivBHChunks(chunks,law,m,wantTangent);end
end
