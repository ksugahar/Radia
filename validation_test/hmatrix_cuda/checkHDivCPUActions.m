function result=checkHDivCPUActions(meshPath,probePath,options,assemblyMode)
% Independent CPU fixture is validation-only, never GPU setup input.
if nargin<4,assemblyMode='hacapk-cached';end
assert(ismember(assemblyMode,{'hacapk-cached','hacapk-device-dense'}));
assert(radia_cuda_hmatrix('live_count')==0,'Use an isolated validation session');
[probes,~]=radia.gpu.readHDivH5(probePath);
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0
    config=radia.gpu.curvedHDivQuadrature(2);config.geometry_mode='curved';
else
    config=radia.gpu.affineHDivQuadrature(2);config.geometry_mode='affine';
end
config.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
    config.(name{1})=setup.(name{1});
end
extra={};if isfield(options,'cache_bytes'),extra={options.cache_bytes};end
command='create_hdiv_hacapk_cached';if strcmp(assemblyMode,'hacapk-device-dense'),command='create_hdiv_hacapk_device_dense';end
h=radia_cuda_hmatrix(command,config,options.epsilon,options.leaf_size,options.eta,options.max_rank,extra{:});
closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
assert(isequal(size(probes.x),size(probes.expected_action)) && size(probes.x,1)==size(setup.mass,1));
errors=zeros(1,size(probes.x,2));
for k=1:numel(errors)
    radia_cuda_hmatrix('upload',h,probes.x(:,k));radia_cuda_hmatrix('apply',h,1);
    actual=radia_cuda_hmatrix('download',h);expected=probes.expected_action(:,k);
    errors(k)=norm(actual-expected)/max(norm(expected),realmin);
end
assert(all(isfinite(errors)) && all(errors<=1e-8),'Independent CPU B-transpose G B action mismatch');
stats=radia_cuda_hmatrix('stats',h);clear closer
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'count',numel(errors),'seed',20261018,'relative_errors',errors, ...
    'limit',1e-8,'true_relative_residual',[],'residual_applicability','not_applicable_operator', ...
    'live_count',0,'operator_stats',stats,'timing_scope','Excluded from workflow benchmarks');
end
