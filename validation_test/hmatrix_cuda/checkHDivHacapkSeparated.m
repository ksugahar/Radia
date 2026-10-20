function result=checkHDivHacapkSeparated(meshPath,outputPath)
% Two translated public affine bodies; operator compression diagnostic only.
arguments
    meshPath (1,:) char
    outputPath (1,:) char=''
end
assert(radia_cuda_hmatrix('live_count')==0,'Use an isolated numerical session');
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
assert(setup.curved_cells==0 && size(setup.cell_nodes,1)==1,'Requires one affine tetrahedron');
config=radia.gpu.affineHDivQuadrature(2);config.geometry_mode='affine';config.schema_version=1;
nc=size(setup.cell_nodes,1);nf=size(setup.face_nodes,1);n=numel(setup.host);
assert(n<=64,'Small diagnostic budget exceeded');
cells=setup.cell_nodes;faces=setup.face_nodes;
cells(:,1:3:end)=cells(:,1:3:end)+8;faces(:,1:3:end)=faces(:,1:3:end)+8;
config.cell_nodes=[setup.cell_nodes;cells];config.face_nodes=[setup.face_nodes;faces];
offset=max([setup.cell_vertices(:);setup.face_vertices(:)])+1;
config.cell_vertices=[setup.cell_vertices;setup.cell_vertices+offset];
config.face_vertices=[setup.face_vertices;setup.face_vertices+offset];
shift=nc*ones(size(setup.host));shift(setup.kind==1)=nf;
config.host=[setup.host;setup.host+shift];config.kind=[setup.kind;setup.kind];config.expo=[setup.expo;setup.expo];
config.charge_map=speye(2*n);
a=radia_cuda_hmatrix('create_hdiv_elements',config);ca=onCleanup(@()radia_cuda_hmatrix('destroy',a));
% Small leaves exercise actual far leaves; this is not a recommended leaf size.
b=radia_cuda_hmatrix('create_hdiv_hacapk',config,1e-10,8,2.,32);cb=onCleanup(@()radia_cuda_hmatrix('destroy',b));
A=zeros(2*n);B=zeros(2*n);
for column=1:2*n
    x=zeros(2*n,1);x(column)=1;
    radia_cuda_hmatrix('upload',a,x);radia_cuda_hmatrix('apply',a,1);A(:,column)=radia_cuda_hmatrix('download',a);
    radia_cuda_hmatrix('upload',b,x);radia_cuda_hmatrix('apply',b,1);B(:,column)=radia_cuda_hmatrix('download',b);
end
error=norm(A-B,'fro')/norm(A,'fro');assert(error<=1e-8,'Separated compressed operator mismatch');
stats=radia_cuda_hmatrix('stats',b);assert(stats.assembly_lowrank_count>0,'No low-rank leaf was exercised');
assert(stats.stored_coefficient_count>0 && stats.same_leaf_dense_coefficient_count>0);
assert(stats.max_leaf_rank>0 && stats.leaf_rank_sum>=stats.max_leaf_rank);
assert(stats.storage_saving_lowrank_leaf_count<=stats.assembly_lowrank_count);
if ~isempty(outputPath)
    radia.gpu.writeHDivH5(outputPath,struct('direct_operator',A,'compressed_operator',B),struct('units','SI','labels',{{'magnetic_body','exterior'}}));
end
clear ca cb
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'charges',2*n,'relative_frobenius_error',error,'live_count',0, ...
    'leaf_size',8,'epsilon',1e-10,'eta',2.,'max_rank',32,'assembly_stats',stats, ...
    'scope','Two translated affine public tetrahedra, GPU operator compression only; not independent CPU, full solve or performance acceptance', ...
    'operation_kind','operator','true_relative_residual',[],'residual_applicability','not_applicable_operator');
end
