function result=checkHDivSparseMass(meshPath)
% Exact storage parity: preserve element accumulation, numbering and quadrature.
a=radia_cuda_projection('assemble_hdiv_elements',meshPath,2,[1 .2 .1]);
b=radia_cuda_projection('assemble_hdiv_elements_sparse_mass',meshPath,2,[1 .2 .1]);
assert(~issparse(a.mass) && issparse(b.mass));
assert(isequal(a.mass,full(b.mass)),'Sparse mass changed FE arithmetic');
for field=fieldnames(a)'
    if ~strcmp(field{1},'mass'),assert(isequal(a.(field{1}),b.(field{1})),'Non-mass staging field changed');end
end
c=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
assert(issparse(c.mass) && issparse(c.charge_map));
assert(isequal(full(c.mass),a.mass) && isequal(full(c.charge_map),a.charge_map));
x=radia.gpu.solveHDivElements(meshPath,2,[1 .2 .1],2.5,1e-11,1000,'host-cache','mass-cholesky','dense');
y=radia.gpu.solveHDivElements(meshPath,2,[1 .2 .1],2.5,1e-11,1000,'host-cache','mass-cholesky','sparse-mass');
z=radia.gpu.solveHDivElements(meshPath,2,[1 .2 .1],2.5,1e-11,1000,'host-cache','mass-cholesky','sparse');
assert(isequal(z.coefficients,x.coefficients));
assert(isequal(x.coefficients,y.coefficients) && isequal(x.mean_magnetization,y.mean_magnetization));
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'exact_mass_equal',true,'exact_staging_equal',true,'exact_solution_equal',true, ...
    'dofs',numel(a.rhs),'mass_nnz',nnz(b.mass),'dense_mass_array_bytes',8*numel(a.mass), ...
    'sparse_mass_payload_bytes',16*nnz(b.mass)+8*(size(b.mass,2)+1),'memory_scope','Numerical array payload; excludes headers, allocator and runtime peaks','live_count',0);
end
