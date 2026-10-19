function result=measureMaterialCSC(mesh,cpu,previousDirectory,candidateDirectory)
% Compare complete nonlinear APIs; the only changed kernel is host CSC staging.
[reference,~]=radia.gpu.readHDivH5(cpu);device=gpuDevice;before=device.AvailableMemory;
directories={previousDirectory,candidateDirectory};samples={};first={};
opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',4000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true,'Preconditioner','mass-polynomial'};
for index=1:2,first{index}=run(index);end
for cycle=1:2
    for index=[1 2 2 1]
        row=run(index);row.cycle=cycle;samples{end+1}=row;
    end
end
samples=[samples{:}];old=samples([samples.variant]==1);candidate=samples([samples.variant]==2);
result=struct('accepted',true,'first_calls',{first},'samples',samples,'previous_whole_median_seconds',median([old.whole_seconds]),'candidate_whole_median_seconds',median([candidate.whole_seconds]),'previous_fe_median_seconds',median([old.fe_and_law_seconds]),'candidate_fe_median_seconds',median([candidate.fe_and_law_seconds]),'previous_over_candidate',median([old.whole_seconds])/median([candidate.whole_seconds]),'available_device_memory_before_sample',before,'available_device_memory_after_sample',device.AvailableMemory,'live_count',radia_cuda_hmatrix('live_count'),'scope','Separate first call per MEX; warmed two-cycle ABBA whole MATLAB API including FE, geometry, transfers, nonlinear solve, recovery and destruction; path selection, clear/reload, Engine, HDF5 and validation excluded; not process cold','memory_limitations','Controlled array bound and available-device-memory samples; host, CUDA library and global peaks unmeasured');
assert(result.live_count==0);
    function row=run(index)
        assert(radia_cuda_hmatrix('live_count')==0);
        for d=directories,if contains(path,d{1}),rmpath(d{1});end,end
        clear radia_cuda_projection;addpath(directories{index},'-begin');
        assert(strcmp(which('radia_cuda_projection'),fullfile(directories{index},'radia_cuda_projection.mexw64')));
        s=radia.gpu.solveNonlinearHDivElements(mesh,2,[0 0 1e5],reference.table,opts{:});
        err=norm(s.coefficients-reference.coefficients(:))/norm(reference.coefficients(:));
        assert(err<=1e-8 && s.nonlinear_relative_residual<=1e-12 && s.live_count==0);
        row=struct('variant',index,'whole_seconds',s.whole_seconds,'prepare_seconds',s.prepare_seconds,'fe_and_law_seconds',s.prepare_phase_seconds.fe_and_law,'solve_seconds',s.solve_seconds,'cpu_coefficient_relative_error',err,'nonlinear_relative_residual',s.nonlinear_relative_residual,'nonlinear_iterations',s.nonlinear_iterations,'inner_iterations',sum([s.history.inner_iterations]),'combined_array_bound_bytes',s.combined_array_bound_bytes,'maximum_solver_owned_device_bytes',s.maximum_solver_owned_device_bytes);
    end
end
