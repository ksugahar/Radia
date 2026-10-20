function result=measureStrongNonlinear(mesh,cpu)
% One accepted configuration: absolute whole timing, not comparative speedup.
[reference,~]=radia.gpu.readHDivH5(cpu);field=reference.applied_field(:)';device=gpuDevice;
opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true,'Preconditioner','mass-chebyshev','FinalResidualOnly',true};
before=device.AvailableMemory;first=record(radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:}));samples={};
for k=1:3,samples{end+1}=record(radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:}));end
samples=[samples{:}];timer=tic;obj=radia.gpu.PreparedHDivNonlinear(mesh,2,field,reference.table,opts{:});owner=onCleanup(@()delete(obj));info=obj.stats();runs={};
for k=1:3,runs{end+1}=record(obj.solve());end
runs=[runs{:}];release=tic;clear owner obj;wait(device);releaseSeconds=toc(release);total=toc(timer);
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'applied_field_A_per_m',field,'first_protocol_sample',first,'whole_samples',samples,'whole_median_seconds',median([samples.whole_seconds]),'prepared_info',info,'prepared_runs',runs,'prepared_solve_median_seconds',median([runs.whole_seconds]),'prepared_construction_three_solves_release_seconds',total,'release_seconds',releaseSeconds,'available_device_memory_before_sample',before,'available_device_memory_after_release_sample',device.AvailableMemory,'live_count',0,'scope','First protocol sample after numerical checks; three warmed complete APIs and three zero-start prepared solves; whole includes FE, H-matrix, transfers, nonlinear solution, recovery and destruction; prepared solve excludes preparation/destruction; Engine/HDF5/validation excluded; no process cold, ABBA or comparative speedup claim','physical_scope','Public synthetic BH with declared finite barrier continuation above the final magnetization knot; not measured steel','memory_limitations','Native controlled bounds and available-memory samples only; library, host and global peaks unmeasured');
    function row=record(s)
        err=norm(s.coefficients-reference.coefficients(:))/norm(reference.coefficients(:));assert(err<=1e-8 && s.nonlinear_relative_residual<=1e-12);
        row=struct('whole_seconds',s.whole_seconds,'prepare_seconds',s.prepare_seconds,'solve_seconds',s.whole_seconds,'nonlinear_updates',s.nonlinear_iterations,'initial_linear_iterations',s.warmstart_iterations,'inner_linear_iterations',sum([s.history.inner_iterations]),'line_search_backtracks',sum([s.history.backtracks]),'history',s.history,'outer_tangent_evaluations',s.outer_tangent_evaluations,'residual_only_outer_evaluations',s.residual_only_outer_evaluations,'cpu_coefficient_relative_error',err,'nonlinear_relative_residual',s.nonlinear_relative_residual,'maximum_material_magnetization_A_per_m',s.maximum_material_magnetization,'barrier_cap_A_per_m',s.material_barrier_cap,'combined_array_bound_bytes',s.combined_array_bound_bytes,'maximum_solver_owned_device_bytes',s.maximum_solver_owned_device_bytes);
        if isfield(s,'solve_seconds'),row.solve_seconds=s.solve_seconds;end
    end
end
