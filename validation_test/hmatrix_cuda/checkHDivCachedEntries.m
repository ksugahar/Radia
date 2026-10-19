function result=checkHDivCachedEntries(meshPath,outputPath,labels)
% Reuse/eviction parity versus the previously accepted mapped-moment dense path.
arguments
    meshPath (1,:) char
    outputPath (1,:) char=''
    labels cell={'magnetic_body','exterior'}
end
assert(radia_cuda_hmatrix('live_count')==0,'Use an isolated numerical session');
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0
    config=radia.gpu.curvedHDivQuadrature(2);config.geometry_mode='curved';
else
    config=radia.gpu.affineHDivQuadrature(2);config.geometry_mode='affine';
end
config.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo'}
    config.(name{1})=setup.(name{1});
end
n=numel(config.host);assert(n>=4 && n<=256,'Small complete-operator budget required');config.charge_map=speye(n);
reference=config;reference.host_vector_cache=true;reference.host_moment_cache=true;
h=radia_cuda_hmatrix('create_hdiv_elements',reference);closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
A=zeros(n);
for column=1:n
    x=zeros(n,1);x(column)=1;
    radia_cuda_hmatrix('upload',h,x);radia_cuda_hmatrix('apply',h,1);A(:,column)=radia_cuda_hmatrix('download',h);
end
[rows,columns]=ndgrid(1:n,1:n);
[values,fullStats]=radia_cuda_hmatrix('curved_entry_batch_cached',config,rows(:),columns(:));B=reshape(values,n,n);
points=size(setup.cell_nodes,1)*size(config.tet_rule,1)+size(setup.face_nodes,1)*size(config.tri_rule,1);
budget=2*floor(n/3)*points*8;assert(budget<=256*2^20);
[values,tileStats]=radia_cuda_hmatrix('curved_entry_batch_cached',config,rows(:),columns(:),budget);C=reshape(values,n,n);
[wideValues,wideStats]=radia_cuda_hmatrix('curved_entry_batch_cached',config,rows(:),columns(:),2*2^30);
wideError=norm(A-reshape(wideValues,n,n),'fro')/norm(A,'fro');
assert(wideError<=1e-12 && wideStats.cache_tile_loads==1,'Explicit larger budget changed the operator');
fullError=norm(A-B,'fro')/norm(A,'fro');tileError=norm(A-C,'fro')/norm(A,'fro');
assert(max(fullError,tileError)<=1e-12,'Indexed potential cache changed the operator');
assert(fullStats.cache_tile_loads==1 && fullStats.cache_tile_hits>0);
assert(tileStats.cache_tile_loads>2,'Forced three-or-more-tile cache did not evict');
selected=[1;n;n;1];other=[1;1;n;n];direct=radia_cuda_hmatrix('curved_entry_batch',config,selected,other);
assert(norm(direct-B(sub2ind([n n],selected,other)))<=1e-12*norm(A,'fro'));
for badBudget={0,1,2*2^30+1}
    rejected=false;try,unexpected=radia_cuda_hmatrix('curved_entry_batch_cached',config,1,1,badBudget{1});catch,rejected=true;end %#ok<NASGU>
    assert(rejected,'Invalid cache budget was accepted');
end
empty=radia_cuda_hmatrix('curved_entry_batch_cached',config,zeros(0,1),zeros(0,1));assert(isempty(empty));
% Reject extended-capacity misuse before constructing a device oracle.
for limits={[4609 64 64],[4097 32 64],[4097 64 257]}
    limits=limits{1};pick=mod((0:limits(1)-1)',n)+1;oversized=config;
    for name={'host','kind','expo'}
        data=config.(name{1});oversized.(name{1})=data(pick,:);
    end
    oversized.charge_map=sparse((1:limits(1))',pick,ones(limits(1),1),limits(1),n);
    rejected=false;
    try
        unexpected=radia_cuda_hmatrix('create_hdiv_hacapk_cached',oversized,1e-10,limits(2),2,limits(3));
        radia_cuda_hmatrix('destroy',unexpected);
    catch failure
        rejected=contains(failure.message,'Unsupported curved HACApK budget');
    end
    assert(rejected,'Extended HACApK dimension/tree/rank guard failed');
end
assert(radia_cuda_hmatrix('live_count')==1);
if ~isempty(outputPath)
    radia.gpu.writeHDivH5(outputPath,struct('reference',A,'cached',B,'evicted',C),struct('units','SI','labels',{labels}));
end
clear closer
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'charges',n,'cached_relative_error',fullError,'evicted_relative_error',tileError, ...
    'wide_relative_error',wideError,'wide_cache',wideStats,'default_cache',fullStats,'forced_eviction_cache',tileStats,'live_count',0,'operation_kind','operator', ...
    'true_relative_residual',[],'residual_applicability','not_applicable_operator', ...
    'scope','GPU full operator/cache eviction parity only; independent CPU and compressed whole solve remain separate');
end
