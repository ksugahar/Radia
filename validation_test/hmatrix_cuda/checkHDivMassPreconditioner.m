function result=checkHDivMassPreconditioner(meshPath)
% GPU factor failure must release work buffers and leave the resident operator usable.
setup=radia_cuda_projection('assemble_hdiv_elements',meshPath,2,[1 .2 .1]);
config=radia.gpu.curvedHDivQuadrature(2);config.geometry_mode='curved';config.schema_version=1;config.host_vector_cache=true;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
    config.(name{1})=setup.(name{1});
end
h=radia_cuda_hmatrix('create_hdiv_elements',config);closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
diagonal=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,2.5,1e-11,1000,'mass-diagonal');
cholesky=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,2.5,1e-11,1000,'mass-cholesky');
assert(norm(diagonal.x-cholesky.x)/norm(cholesky.x)<1e-8);
assert(max([diagonal.true_relative_residual cholesky.true_relative_residual])<1e-10);
bad=eye(numel(setup.rhs));bad(1,2)=2;bad(2,1)=2;rejected=false;
try
    unexpected=radia_cuda_hmatrix('solve_linear_material',h,bad,setup.rhs,2.5,1e-11,1000,'mass-cholesky');
catch exception
    rejected=contains(exception.message,'non-SPD');
end
assert(rejected,'Expected explicit non-SPD factor rejection');assert(radia_cuda_hmatrix('live_count')==1);
repeat=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,2.5,1e-11,1000,'mass-cholesky');
assert(isequal(repeat.x,cholesky.x),'Failed factor changed repeat solution');
zero=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,zeros(size(setup.rhs)),2.5,1e-11,1000,'mass-cholesky');
assert(zero.iterations==0 && zero.true_relative_residual==0 && all(zero.x==0));
result=struct('accepted',true,'diagonal_iterations',diagonal.iterations,'mass_cholesky_iterations',cholesky.iterations,'non_spd_rejected',rejected,'repeat_unchanged',true,'zero_rhs_passed',true);
clear closer
assert(radia_cuda_hmatrix('live_count')==0);result.live_count=0;
end
