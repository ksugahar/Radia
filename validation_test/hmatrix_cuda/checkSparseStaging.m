function result=checkSparseStaging(mesh,previousDirectory,candidateDirectory)
% Compare FE mass/charge staging after allocation-admission changes.
commands={'assemble_hdiv_elements_sparse','assemble_hdiv_elements_sparse_mass'};
old=cell(size(commands));
addpath(previousDirectory,'-begin');clear radia_cuda_projection;
assert(strcmp(which('radia_cuda_projection'),fullfile(previousDirectory,'radia_cuda_projection.mexw64')));
for k=1:numel(commands),old{k}=radia_cuda_projection(commands{k},mesh,2,[0 0 1e5]);end
assert(radia_cuda_projection('live_count')==0);
rmpath(previousDirectory);clear radia_cuda_projection;addpath(candidateDirectory,'-begin');
assert(strcmp(which('radia_cuda_projection'),fullfile(candidateDirectory,'radia_cuda_projection.mexw64')));
for k=1:numel(commands)
    actual=radia_cuda_projection(commands{k},mesh,2,[0 0 1e5]);
    assert(isequal(old{k},actual),'FE staging fields or values changed');
    assert(issparse(actual.mass));assert(issparse(actual.charge_map)==(k==1));
end
assert(isequal(old{1}.mass,old{2}.mass));
assert(isequal(old{1}.charge_map,old{2}.charge_map));
rejected=0;
for order=[0 4]
    try,unused=radia_cuda_projection(commands{1},mesh,order,[0 0 1e5]);catch,rejected=rejected+1;end
end
assert(rejected==2);again=radia_cuda_projection(commands{1},mesh,2,[0 0 1e5]);assert(isequal(old{1},again));
assert(radia_cuda_projection('live_count')==0);
result=struct('accepted',true,'dofs',numel(again.rhs),'mass_nonzeros',nnz(again.mass),'charge_nonzeros',nnz(again.charge_map),'all_fields_exact',true,'invalid_orders_rejected',rejected,'live_count',0,'scope','Exact FE staging regression; no new solver or performance acceptance');
clear old actual again;
rmpath(candidateDirectory);clear radia_cuda_projection;
end
