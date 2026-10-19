function r=checkCurvedAllocationFailure(meshPath,referencePath,mode,kind,target,site)
% Public curved P1 fixture. Native setup and CPU reference are independently checked.
if nargin<6,site='';end
persistent cachedMesh cachedReference cfg setup reference setupErrors
if isempty(cachedMesh)||~strcmp(cachedMesh,meshPath)||~strcmp(cachedReference,referencePath)
    reference=load(referencePath);
    setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,1,[1 .2 .1]);
    assert(setup.curved_cells>0,'Expected the frozen curved mesh');
    setupErrors=[norm(setup.mass-full(reference.mass),'fro')/norm(full(reference.mass),'fro'), ...
        norm(setup.charge_map-full(reference.B),'fro')/norm(full(reference.B),'fro'), ...
        norm(setup.rhs-reference.rhs)/norm(reference.rhs)];
    if max(setupErrors)>=1e-11
        r=struct('accepted',false,'mode',mode,'kind',kind,'target',target,'setup_relative_errors',setupErrors,'error','Independent CPU FE setup rejected');return;
    end
    cfg=radia.gpu.curvedHDivQuadrature(1);cfg.geometry_mode='curved';cfg.schema_version=1;
    for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
        cfg.(name{1})=setup.(name{1});
    end
    cachedMesh=meshPath;cachedReference=referencePath;
end
switch mode
    case 'direct',command='create_hdiv_hacapk';extra={};
    case 'cached',command='create_hdiv_hacapk_cached';extra={32*2^20};
    case 'device-dense',command='create_hdiv_hacapk_device_dense';extra={32*2^20};
    case {'resident','resident-cache','resident-moments','resident-profile'}
        command='create_curved_hdiv';extra={};
    otherwise,error('Unsupported curved diagnostic mode');
end
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_resource_failure',0);
radia_cuda_hmatrix('diagnostic_publication_failure',0);
radia_cuda_hmatrix('diagnostic_trapped_frame_exception',0);
switch kind
    case 'c',radia_cuda_hmatrix('diagnostic_c_allocation_failure',target);
    case {'cpp','cpp-baseline'}
        if isempty(site),radia_cuda_hmatrix('diagnostic_host_allocation_failure',target);
        else,radia_cuda_hmatrix('diagnostic_host_allocation_failure',target,site);end
    case 'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_failure',target);
    case 'publication',radia_cuda_hmatrix('diagnostic_publication_failure',target);
    case 'resource',radia_cuda_hmatrix('diagnostic_resource_failure',target);
    case 'baseline'
    otherwise,error('Unsupported curved diagnostic failure class');
end
h=[];failed=false;identifier='';message='';numerical=struct();createCounts=struct();
try
    if startsWith(mode,'resident')
        localConfig=cfg;
        localConfig.host_vector_cache=~strcmp(mode,'resident');
        localConfig.host_moment_cache=any(strcmp(mode,{'resident-moments','resident-profile'}));
        localConfig.profile_assembly=strcmp(mode,'resident-profile');
        h=radia_cuda_hmatrix(command,localConfig);
    else
        h=radia_cuda_hmatrix(command,cfg,1e-12,64,2.,64,extra{:});
    end
    % Enumerate this create route before the baseline solve adds its own
    % allocations. Failure injection below deliberately tests create only.
    createC=radia_cuda_hmatrix('diagnostic_c_allocation_report');
    createCpp=radia_cuda_hmatrix('diagnostic_host_allocation_report');
    createGpu=radia_cuda_hmatrix('diagnostic_device_allocation_report');
    createResource=radia_cuda_hmatrix('diagnostic_resource_report');
    createCounts=struct('c',createC.allocations,'cpp',createCpp.attempts, ...
        'gpu',createGpu.attempts,'resource',createResource.attempts);
    if strcmp(kind,'baseline')
        s=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,2.5,1e-11,1000);
        numerical=struct('cpu_error',norm(s.x-reference.expected)/norm(reference.expected), ...
            'physical_error',norm(setup.mean_operator*s.x-reference.mean_operator*reference.expected)/norm(reference.mean_operator*reference.expected), ...
            'independent_true_residual',norm(reference.rhs-reference.system*s.x)/norm(reference.rhs), ...
            'solver_true_residual',s.true_relative_residual,'iterations',s.iterations);
        numerical.accepted=max([numerical.cpu_error numerical.physical_error])<1e-8 && numerical.independent_true_residual<1e-10 && numerical.solver_true_residual<=1e-11;
        stats=radia_cuda_hmatrix('stats',h);
        if strcmp(mode,'device-dense'),numerical.accepted=numerical.accepted && stats.assembly_deferred_dense_leaves>0;end
    end
catch e
    failed=true;identifier=e.identifier;message=e.message;
end
if ~isempty(h),radia_cuda_hmatrix('destroy',h);end
r=struct('mode',mode,'kind',kind,'target',target,'allocation_site',site,'failed',failed,'identifier',identifier, ...
    'message',message,'numerical',numerical,'create_counts',createCounts,'setup_relative_errors',setupErrors, ...
    'c',radia_cuda_hmatrix('diagnostic_c_allocation_report'), ...
    'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_report'), ...
    'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_report'), ...
    'resources',radia_cuda_hmatrix('diagnostic_resource_report'), ...
    'publication',radia_cuda_hmatrix('diagnostic_publication_report'), ...
    'live_count',radia_cuda_hmatrix('live_count'));
clean=r.live_count==0 && r.c.live_bytes==0 && r.cpp.live_bytes==0 && r.gpu.live_bytes==0 && r.gpu.free_errors==0;
clean=clean && r.resources.live_resources==0 && r.resources.free_errors==0 && r.publication.owned_output_arrays==0;
if strcmp(kind,'baseline'),expected=~failed && isfield(numerical,'accepted') && numerical.accepted;
elseif strcmp(kind,'cpp-baseline'),expected=~failed && target==0;
else,expected=failed && strcmp(identifier,'radia:cuda:hmatrix');end
if strcmp(kind,'c'),expected=expected && r.c.allocations==target && r.c.failure_code==1;end
if strcmp(kind,'cpp'),expected=expected && r.cpp.attempts==target;end
if ~isempty(site),expected=expected && strcmp(r.cpp.filter,site);end
if strcmp(kind,'gpu'),expected=expected && r.gpu.attempts==target;end
if strcmp(kind,'resource'),expected=expected && r.resources.attempts==target;end
if strcmp(kind,'publication'),expected=expected && r.publication.reached_stage==target;end
r.accepted=clean && expected;
if ~clean,return;end
% Release diagnostic selection; the cached MATLAB input is small fixture data only.
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_resource_failure',0);
radia_cuda_hmatrix('diagnostic_publication_failure',0);
end
