function result=solveCurvedHDiv(meshPath,order,appliedField,invChi,tolerance,maxIterations)
% Experimental bounded P2 tetrahedral HDiv material solve, no Python setup.
arguments
    meshPath (1,:) char
    order (1,1) double {mustBeMember(order,[1 2])}
    appliedField (1,3) double {mustBeFinite}
    invChi (1,1) double {mustBePositive,mustBeFinite}
    tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
    maxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
end
assert(tolerance<=1e-6,'True residual tolerance must be at most 1e-6');
timer=tic;
setup=radia_cuda_projection('assemble_hdiv_elements',meshPath,order,appliedField);
assert(setup.curved_cells>0,'Curved reference-frame assembly requires a curved mesh');
config=radia.gpu.curvedHDivQuadrature(order);
config.schema_version=1;
config.geometry_mode='curved';
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
    config.(name{1})=setup.(name{1});
end
h=radia_cuda_hmatrix('create_curved_hdiv',config);
closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
solution=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,invChi,tolerance,maxIterations);
stats=radia_cuda_hmatrix('stats',h);
result=struct('coefficients',solution.x,'mean_magnetization',setup.mean_operator*solution.x, ...
    'true_relative_residual',solution.true_relative_residual,'iterations',solution.iterations, ...
    'operator_owned_device_bytes',stats.owned_device_bytes,'solver_owned_device_bytes',solution.solver_owned_device_bytes, ...
    'setup_peak_owned_device_bytes',stats.setup_peak_owned_device_bytes, ...
    'combined_peak_owned_device_bytes',max(stats.setup_peak_owned_device_bytes,stats.owned_device_bytes+solution.solver_owned_device_bytes), ...
    'scope','Bounded P2 tetrahedral reference-monomial dense Gram, linear isotropic material');
clear closer
result.whole_seconds=toc(timer);
result.memory_limitations='Accounted owned buffers; no cuBLAS, host or global allocation peak.';
end
