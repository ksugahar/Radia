function result=prepareCurvedMxOutputFixture(meshPath,cpuPath,outputPath)
% CPU FE operator is independent; charge-entry comparisons remain GPU-internal.
reference=load(cpuPath);
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,1,[1 .2 .1]);
assert(setup.curved_cells>0);
setupErrors=[norm(setup.mass-full(reference.mass),'fro')/norm(full(reference.mass),'fro'), ...
    norm(setup.charge_map-full(reference.B),'fro')/norm(full(reference.B),'fro'), ...
    norm(setup.rhs-reference.rhs)/norm(reference.rhs)];
assert(max(setupErrors)<=1e-11);
config=radia.gpu.curvedHDivQuadrature(1);config.geometry_mode='curved';config.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
    config.(name{1})=setup.(name{1});
end
n=numel(config.host);assert(n<=128,'Small complete charge operator fixture only');
[row,column]=ndgrid(1:n,1:n);
values=radia_cuda_hmatrix('curved_entry_batch',config,row(:),column(:));
gram=reshape(values,n,n);
actual=2.5*full(reference.mass)+full(reference.B)'*gram*full(reference.B);
cpuError=norm(actual-reference.system,'fro')/norm(reference.system,'fro');
assert(cpuError<=1e-10,'Independent CPU mapped operator rejected');
rows=[1;n;1;n];columns=[1;1;n;n];
curved_output_probe=struct('config',config,'rows',rows,'columns',columns, ...
    'expected',gram(sub2ind([n n],rows,columns)), ...
    'cache_budget_bytes',32*2^20,'cpu_operator_relative_error',cpuError, ...
    'setup_relative_errors',setupErrors);
save(outputPath,'curved_output_probe','-v7');
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'charges',n,'cpu_operator_relative_error',cpuError, ...
    'setup_relative_errors',setupErrors,'scope', ...
    'Independent CPU mapped material operator; per-entry reference is direct GPU, not independent charge energy');
end
