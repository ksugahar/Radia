function result=checkMaterialCSC(mesh,previousDirectory,candidateDirectory)
% Exact old/new staging comparison; native regions timed separately, no speed claim.
addpath(previousDirectory,'-begin');clear radia_cuda_projection;
assert(strcmp(which('radia_cuda_projection'),fullfile(previousDirectory,'radia_cuda_projection.mexw64')));
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',mesh,2,[0 0 1e5]);
ne=size(setup.cell_nodes,1);old={};oldSeconds=0;
for first=0:128:ne-1
    timer=tic;old{end+1}=radia_cuda_projection('hdiv_material_quadrature_chunk',mesh,2,first,min(128,ne-first));oldSeconds=oldSeconds+toc(timer);
end
rmpath(previousDirectory);clear radia_cuda_projection;addpath(candidateDirectory,'-begin');
assert(strcmp(which('radia_cuda_projection'),fullfile(candidateDirectory,'radia_cuda_projection.mexw64')));
newSeconds=0;nonzeros=0;fields={'weights','points','element','global_point_indices','integration_order','dofs','total_points','total_elements','first_element','element_count','potential_shape_entries','tangent_nonzero_bound'};
for k=1:numel(old)
    first=(k-1)*128;timer=tic;q=radia_cuda_projection('hdiv_material_quadrature_chunk',mesh,2,first,min(128,ne-first));newSeconds=newSeconds+toc(timer);
    a=old{k};assert(isequal(a.evaluation,q.evaluation),'Evaluation entries differ');
    [ai,aj,av]=find(a.evaluation);[bi,bj,bv]=find(q.evaluation);assert(isequal(ai,bi)&&isequal(aj,bj)&&isequal(av,bv),'Sparse structure differs');
    for f=fields,assert(isequal(a.(f{1}),q.(f{1})),f{1});end
    nonzeros=nonzeros+nnz(q.evaluation);old{k}=[];
end
rejected=0;
for args={[-1 1],[0 0],[ne 1],[0 ne+1]}
    v=args{1};try,q=radia_cuda_projection('hdiv_material_quadrature_chunk',mesh,2,v(1),v(2));catch,rejected=rejected+1;end
end
assert(rejected==4);
result=struct('accepted',true,'elements',ne,'dofs',numel(setup.rhs),'chunks',numel(old),'nonzeros',nonzeros,'evaluation_bitwise_equal',true,'metadata_equal',true,'invalid_intervals_rejected',rejected,'previous_staging_seconds',oldSeconds,'candidate_staging_seconds',newSeconds,'scope','Ordered single-pass host staging diagnostic; not ABBA or GPU whole performance; temporary arrays released on return');
end
