function result=solveHDivRT0(meshPath,appliedField,invChi,tolerance,maxIterations,assembly)
% Experimental affine tetra RT0 material solve, fully MATLAB/native setup.
arguments
    meshPath (1,:) char
    appliedField (1,3) double
    invChi (1,1) double {mustBePositive,mustBeFinite}
    tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
    maxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
    assembly (1,:) char {mustBeMember(assembly,{'dense','hacapk'})}='dense'
end
assert(tolerance<=1e-6,'True residual tolerance must be at most 1e-6');
assert(all(isfinite(appliedField)),'Applied field must be finite');
timer=tic;
setup=radia_cuda_projection('assemble_hdiv_rt0',meshPath,appliedField);
config=struct('schema_version',1,'operator_kind','demag','assembly','analytic_constant_block', ...
    'cells',setup.cells,'faces',setup.faces,'charge_map',setup.charge_map);
if strcmp(assembly,'hacapk')
    config=rmfield(config,'assembly');
    config.epsilon=1e-12;config.leaf_size=64;config.eta=2;config.max_rank=32;
    h=radia_cuda_hmatrix('create_analytic_hacapk',config);
else
    h=radia_cuda_hmatrix('create_analytic_constant_block',config);
end
closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
solution=radia_cuda_hmatrix('solve_linear_material',h,setup.mass,setup.rhs,invChi,tolerance,maxIterations);
stats=radia_cuda_hmatrix('stats',h);
result=struct('coefficients',solution.x,'mean_magnetization',setup.mean_operator*solution.x, ...
    'true_relative_residual',solution.true_relative_residual,'iterations',solution.iterations, ...
    'preconditioner',solution.preconditioner,'operator_owned_device_bytes',stats.owned_device_bytes, ...
    'solver_owned_device_bytes',solution.solver_owned_device_bytes, ...
    'setup_peak_owned_device_bytes',stats.setup_peak_owned_device_bytes, ...
    'combined_peak_owned_device_bytes',max(stats.setup_peak_owned_device_bytes,stats.owned_device_bytes+solution.solver_owned_device_bytes), ...
    'assembly',assembly,'scope','Bounded affine-tetra RT0, analytic constant-charge Gram, linear isotropic material');
clear closer
result.whole_seconds=toc(timer);
result.memory_limitations='Owned operator and solver buffers; no process/global peak. Setup scratch and solver scratch have different lifetimes.';
end
