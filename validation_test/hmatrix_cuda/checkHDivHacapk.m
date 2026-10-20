function result=checkHDivHacapk(meshPath,leafSize)
% GPU route parity only; independent CPU and scaling gates remain separate.
arguments
    meshPath (1,:) char
    leafSize (1,1) double {mustBeInteger,mustBePositive}=64
end
assert(radia_cuda_hmatrix('live_count')==0,'Use an isolated numerical session');
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
assert(numel(config.host)<=256,'Bounded diagnostic budget exceeded');
direct=radia_cuda_hmatrix('create_hdiv_elements',config);cd=onCleanup(@()radia_cuda_hmatrix('destroy',direct));
compressed=radia_cuda_hmatrix('create_hdiv_hacapk',config,1e-10,leafSize,2.,64);cc=onCleanup(@()radia_cuda_hmatrix('destroy',compressed));
n=size(setup.mass,1);relative=0;
for k=1:3
    x=sin((1:n)'*(.13+.17*k));
    radia_cuda_hmatrix('upload',direct,x);radia_cuda_hmatrix('apply',direct,1);a=radia_cuda_hmatrix('download',direct);
    radia_cuda_hmatrix('upload',compressed,x);radia_cuda_hmatrix('apply',compressed,1);b=radia_cuda_hmatrix('download',compressed);
    relative=max(relative,norm(a-b)/max(norm(a),realmin));
end
assert(relative<=1e-8,'Compressed B-transpose G B action mismatch');
a=radia_cuda_hmatrix('solve_linear_material',direct,setup.mass,setup.rhs,2.5,1e-11,1000,'mass-cholesky');
b=radia_cuda_hmatrix('solve_linear_material',compressed,setup.mass,setup.rhs,2.5,1e-11,1000,'mass-cholesky');
solutionError=norm(a.x-b.x)/norm(a.x);assert(solutionError<=1e-8,'Compressed material solution mismatch');
radia_cuda_hmatrix('upload',direct,b.x);radia_cuda_hmatrix('apply',direct,1);
originalAction=radia_cuda_hmatrix('download',direct)+2.5*(setup.mass*b.x);
originalResidual=norm(originalAction-setup.rhs)/norm(setup.rhs);
assert(originalResidual<=1e-11,'Original direct-operator equation residual failed');
stats=radia_cuda_hmatrix('stats',compressed);
clear cd cc
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'action_relative_error',relative,'solution_relative_error',solutionError, ...
    'true_relative_residual',b.true_relative_residual,'original_direct_true_relative_residual',originalResidual, ...
    'leaf_size',leafSize,'assembly_stats',stats,'live_count',0, ...
    'scope','Bounded direct/compressed GPU parity only; independent CPU original-operator residual, low-rank compression and performance acceptance still required');
end
