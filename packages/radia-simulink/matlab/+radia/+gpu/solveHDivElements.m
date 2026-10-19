function result=solveHDivElements(meshPath,order,appliedField,invChi,tolerance,maxIterations,assemblyMode,preconditioner,stagingMode,profileAssembly,hacapkOptions)
% Experimental bounded P2 tetrahedral HDiv material solve, no Python setup.
arguments
    meshPath (1,:) char
    order (1,1) double {mustBeMember(order,[1 2])}
    appliedField (1,3) double {mustBeFinite}
    invChi (1,1) double {mustBePositive,mustBeFinite}
    tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
    maxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
    assemblyMode (1,:) char {mustBeMember(assemblyMode,{'direct','host-cache','host-cache-moments','hacapk','hacapk-cached','hacapk-device-dense'})}='direct'
    preconditioner (1,:) char {mustBeMember(preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial'})}='mass-diagonal'
    stagingMode (1,:) char {mustBeMember(stagingMode,{'dense','sparse-mass','sparse'})}='dense'
    profileAssembly (1,1) logical=false
    hacapkOptions (1,1) struct=struct()
end
assert(~profileAssembly || ~strcmp(assemblyMode,'direct'),'Assembly profile requires a cached assembly mode');
if startsWith(assemblyMode,'hacapk')
    assert(~profileAssembly,'HACApK entry profiling is unsupported');
    required={'epsilon';'leaf_size';'eta';'max_rank'};
    provided=fieldnames(hacapkOptions);
    if isfield(hacapkOptions,'cache_bytes')
        assert(ismember(assemblyMode,{'hacapk-cached','hacapk-device-dense'}),'cache_bytes requires cached HACApK');
        validateattributes(hacapkOptions.cache_bytes,{'double'},{'scalar','real','finite','integer','positive','<=',2*2^30});
        provided=setdiff(provided,{'cache_bytes'});
    end
    assert(isequal(sort(provided),sort(required)), ...
        'HACApK requires epsilon, leaf_size, eta, max_rank and optional cached cache_bytes');
else
    assert(isempty(fieldnames(hacapkOptions)),'HACApK options require the explicit hacapk assembly mode');
end
assert(tolerance<=1e-6,'True residual tolerance must be at most 1e-6');
timer=tic;
stageTimer=tic;
if strcmp(stagingMode,'sparse'),command='assemble_hdiv_elements_sparse';elseif strcmp(stagingMode,'sparse-mass'),command='assemble_hdiv_elements_sparse_mass';else,command='assemble_hdiv_elements';end
setup=radia_cuda_projection(command,meshPath,order,appliedField);
feSetupSeconds=toc(stageTimer);stageTimer=tic;
if setup.curved_cells>0
    config=radia.gpu.curvedHDivQuadrature(order);config.geometry_mode='curved';
else
    config=radia.gpu.affineHDivQuadrature(order);config.geometry_mode='affine';
end
config.schema_version=1;
if profileAssembly,config.profile_assembly=true;end
if startsWith(assemblyMode,'host-cache'),config.host_vector_cache=true;end
if strcmp(assemblyMode,'host-cache-moments'),config.host_moment_cache=true;end
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
    config.(name{1})=setup.(name{1});
end
if startsWith(assemblyMode,'hacapk')
    if strcmp(assemblyMode,'hacapk-device-dense'),command='create_hdiv_hacapk_device_dense';elseif strcmp(assemblyMode,'hacapk-cached'),command='create_hdiv_hacapk_cached';else,command='create_hdiv_hacapk';end
    extra={};if isfield(hacapkOptions,'cache_bytes'),extra={hacapkOptions.cache_bytes};end
    h=radia_cuda_hmatrix(command,config,hacapkOptions.epsilon,hacapkOptions.leaf_size,hacapkOptions.eta,hacapkOptions.max_rank,extra{:});
else
    h=radia_cuda_hmatrix('create_hdiv_elements',config);
end
assemblyUploadSeconds=toc(stageTimer);
closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
stageTimer=tic;
solution=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,invChi,tolerance,maxIterations,preconditioner);
solveTransferSeconds=toc(stageTimer);stageTimer=tic;
stats=radia_cuda_hmatrix('stats',h);
assert(~profileAssembly || (isfield(stats,'assembly_profile_enabled') && stats.assembly_profile_enabled),'Loaded native module does not support requested assembly profiling');
result=struct('coefficients',solution.x,'mean_magnetization',setup.mean_operator*solution.x, ...
    'true_relative_residual',solution.true_relative_residual,'iterations',solution.iterations, ...
    'operator_owned_device_bytes',stats.owned_device_bytes,'solver_owned_device_bytes',solution.solver_owned_device_bytes, ...
    'setup_peak_owned_device_bytes',stats.setup_peak_owned_device_bytes, ...
    'combined_peak_owned_device_bytes',max(stats.setup_peak_owned_device_bytes,stats.owned_device_bytes+solution.solver_owned_device_bytes), ...
    'scope','Bounded tetrahedral high-order dense near-rule Gram, linear isotropic material');
recoverySeconds=toc(stageTimer);stageTimer=tic;
clear closer
result.destroy_seconds=toc(stageTimer);
result.whole_seconds=toc(timer);
result.fe_setup_seconds=feSetupSeconds;
result.assembly_upload_seconds=assemblyUploadSeconds;
result.solve_transfer_seconds=solveTransferSeconds;
result.recovery_seconds=recoverySeconds;
result.timing_limitations='Stage wall times include native transfers; transfer-only and kernel-only times are not separated.';
result.assembly_mode=assemblyMode;
result.hacapk_options=hacapkOptions;
if startsWith(assemblyMode,'hacapk'),result.scope='Bounded tetrahedral high-order HACApK Gram, linear isotropic material';end
result.assembly_profile=stats;
result.assembly_profile_limitations='Opt-in CUDA event stream spans include launch gaps and event overhead; diagnostic only, not transfer-inclusive performance or owned-buffer memory. Stage times are zero when profiling is disabled.';
result.staging_mode=stagingMode;
result.preconditioner=solution.preconditioner;
result.memory_limitations='Accounted owned buffers; no cuBLAS, host or global allocation peak.';
end
