function result=checkMaterialContext(mesh,previousDirectory,candidateDirectory)
% Immutable NGSolve mesh/space reuse; identical mapped values and numbering.
addpath(previousDirectory,'-begin');clear radia_cuda_projection;
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',mesh,2,[0 0 1e5]);ne=size(setup.cell_nodes,1);old={};
for first=0:128:ne-1,old{end+1}=radia_cuda_projection('hdiv_material_quadrature_chunk',mesh,2,first,min(128,ne-first));end
assert(radia_cuda_projection('live_count')==0);rmpath(previousDirectory);clear radia_cuda_projection;
addpath(candidateDirectory,'-begin');assert(strcmp(which('radia_cuda_projection'),fullfile(candidateDirectory,'radia_cuda_projection.mexw64')));
h=radia_cuda_projection('create_hdiv_material_context',mesh,2);owner=onCleanup(@()radia_cuda_projection('destroy_hdiv_material_context',h));
assert(radia_cuda_projection('live_count')==1);
for k=1:numel(old)
    first=(k-1)*128;count=min(128,ne-first);
    actual=radia_cuda_projection('hdiv_material_context_chunk',h,first,count);
    direct=radia_cuda_projection('hdiv_material_quadrature_chunk',mesh,2,first,count);
    assert(isequal(old{k},actual)&&isequal(actual,direct),'Material map or metadata changed');old{k}=[];
end
rejected=0;
for interval={[-1 1],[0 0],[ne 1],[0 ne+1],[0 .5]}
    x=interval{1};try,unused=radia_cuda_projection('hdiv_material_context_chunk',h,x(1),x(2));catch,rejected=rejected+1;end
end
try,unused=radia_cuda_projection('points',h);catch,rejected=rejected+1;end
assert(rejected==6 && radia_cuda_projection('live_count')==1);
a=radia_cuda_projection('hdiv_material_context_chunk',h,0,min(128,ne));b=radia_cuda_projection('hdiv_material_context_chunk',h,0,min(128,ne));assert(isequal(a,b));
clear owner;
stale=false;try,unused=radia_cuda_projection('hdiv_material_context_chunk',h,0,1);catch,stale=true;end;assert(stale && radia_cuda_projection('live_count')==0);
clear radia_cuda_projection;
staleReload=false;try,unused=radia_cuda_projection('hdiv_material_context_chunk',h,0,1);catch,staleReload=true;end;assert(staleReload);
fresh=radia_cuda_projection('create_hdiv_material_context',mesh,2);cleanup=onCleanup(@()radia_cuda_projection('destroy_hdiv_material_context',fresh));
assert(fresh~=h);unused=radia_cuda_projection('hdiv_material_context_chunk',fresh,0,1);clear cleanup;
assert(radia_cuda_projection('live_count')==0);
result=struct('accepted',true,'dofs',numel(setup.rhs),'elements',ne,'chunks',numel(old),'all_fields_exact',true,'invalid_and_wrong_type_rejected',rejected,'stale_after_destroy_rejected',stale,'stale_after_reload_rejected',staleReload,'live_count',0,'scope','Native immutable material staging context; no PreparedHDiv integration or performance acceptance');
rmpath(candidateDirectory);clear radia_cuda_projection;
end
