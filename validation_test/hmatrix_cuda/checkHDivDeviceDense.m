function result=checkHDivDeviceDense(meshPath)
% Same-module operator/transaction gate; independent CPU gate is separate.
assert(radia_cuda_hmatrix('live_count')==0);
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0,cfg=radia.gpu.curvedHDivQuadrature(2);cfg.geometry_mode='curved';
else,cfg=radia.gpu.affineHDivQuadrature(2);cfg.geometry_mode='affine';end
cfg.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'},cfg.(name{1})=setup.(name{1});end
n=numel(cfg.host);assert(n<=256,'Bounded full-operator gate');cfg.charge_map=speye(n);
h=radia_cuda_hmatrix('create_hdiv_hacapk_device_dense',cfg,1e-10,64,2.,64,32*2^20);owner=onCleanup(@()radia_cuda_hmatrix('destroy',h));
old=radia_cuda_hmatrix('create_hdiv_hacapk_cached',cfg,1e-10,64,2.,64,32*2^20);oldOwner=onCleanup(@()radia_cuda_hmatrix('destroy',old));
a=zeros(n);b=zeros(n);
for k=1:n
    x=zeros(n,1);x(k)=1;
    radia_cuda_hmatrix('upload',h,x);radia_cuda_hmatrix('apply',h,1);a(:,k)=radia_cuda_hmatrix('download',h);
    radia_cuda_hmatrix('upload',old,x);radia_cuda_hmatrix('apply',old,1);b(:,k)=radia_cuda_hmatrix('download',old);
end
relative=norm(a-b,'fro')/max(norm(b,'fro'),realmin);assert(relative<=1e-11);
stats=radia_cuda_hmatrix('stats',h);
assert(stats.assembly_deferred_dense_leaves>0 && stats.assembly_deferred_dense_coefficients>0);
assert(stats.assembly_coefficient_upload_bytes<stats.stored_coefficient_count*8);
assert(stats.setup_peak_owned_device_bytes>stats.owned_device_bytes && stats.setup_peak_owned_device_bytes<=2*2^30);
baseline=radia_cuda_hmatrix('live_count');
for failure=1:2
    rejected=false;
    try
        if failure==1
            radia_cuda_hmatrix('create_hdiv_hacapk_device_dense',cfg,1e-10,64,2.,64,1);
        else
            radia_cuda_hmatrix('create_hdiv_hacapk_device_dense',cfg,1e-10,64,2.,4097,32*2^20);
        end
    catch,rejected=true;end
    assert(rejected && radia_cuda_hmatrix('live_count')==baseline);
end
clear owner oldOwner
staleRejected=false;try,radia_cuda_hmatrix('apply',h,1);catch,staleRejected=true;end
assert(staleRejected && radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'charges',n,'relative_frobenius_error',relative,'stats',stats, ...
    'failed_construction_count',2,'stale_rejected',staleRejected,'live_count',0, ...
    'scope','Bounded same-module GPU whole-operator and construction/cleanup only; independent CPU acceptance separate');
end
