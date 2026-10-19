function result=measurePolynomialGraph(mesh,cpu)
% All-transfer one-shot ABBA; explicit setup amortization for retained contexts.
[reference,~]=radia.gpu.readHDivH5(cpu);device=gpuDevice;
opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',4000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true};
preconditioners={'mass-polynomial','mass-polynomial-graph'};samples={};initialFree=device.AvailableMemory;
firstOriginal=summarize(radia.gpu.solveNonlinearHDivElements(mesh,2,[0 0 1e5],reference.table,opts{:},'Preconditioner',preconditioners{1}));
firstGraph=summarize(radia.gpu.solveNonlinearHDivElements(mesh,2,[0 0 1e5],reference.table,opts{:},'Preconditioner',preconditioners{2}));
for cycle=1:2
    for mode=[false true true false]
        s=radia.gpu.solveNonlinearHDivElements(mesh,2,[0 0 1e5],reference.table,opts{:},'Preconditioner',preconditioners{1+mode});
        record=summarize(s);record.graph=mode;record.cycle=cycle;samples{end+1}=record;
        assert(s.live_count==0);
    end
end
samples=[samples{:}];series=struct();
for mode=[false true]
    timer=tic;obj=radia.gpu.PreparedHDivNonlinear(mesh,2,[0 0 1e5],reference.table,opts{:},'Preconditioner',preconditioners{1+mode});owner=onCleanup(@()delete(obj));
    info=obj.stats();runs={};
    for k=1:3,runs{end+1}=summarize(obj.solve());end
    runs=[runs{:}];
    retainedFree=device.AvailableMemory;release=tic;clear owner obj;wait(device);releaseSeconds=toc(release);total=toc(timer);
    assert(radia_cuda_hmatrix('live_count')==0);
    name='original';if mode,name='graph';end
    series.(name)=struct('prepare',info,'runs',runs,'release_seconds',releaseSeconds,'construction_three_solves_release_seconds',total,'available_memory_retained_sample',retainedFree);
end
a=[samples(~[samples.graph]).whole_seconds];b=[samples([samples.graph]).whole_seconds];
result=struct('accepted',true,'first_original',firstOriginal,'first_graph',firstGraph,'samples',samples,'one_shot_original_median_seconds',median(a),'one_shot_graph_median_seconds',median(b),'original_to_graph_whole_ratio',median(a)/median(b),'prepared_series',series,'available_device_memory_before_sample',initialFree,'available_device_memory_after_release_sample',device.AvailableMemory,'live_count',radia_cuda_hmatrix('live_count'),'scope','Separate first calls followed by two warmed one-shot ABBA cycles include FE, H-matrix, transfers, Newton, recovery and release; prepared series separately include construction and release; Engine and one-shot validation excluded; prepared-series total includes small result checks; not process cold','memory_limitations','Controlled native/MATLAB bounds and available-device-memory samples only; CUDA graph driver allocations and library/host/global peak unmeasured');
    function r=summarize(s)
        err=norm(s.coefficients-reference.coefficients(:))/norm(reference.coefficients(:));
        assert(err<=1e-8 && s.nonlinear_relative_residual<=1e-12);
        r=struct('whole_seconds',s.whole_seconds,'cpu_coefficient_relative_error',err,'nonlinear_relative_residual',s.nonlinear_relative_residual,'nonlinear_iterations',s.nonlinear_iterations,'warmstart_iterations',s.warmstart_iterations,'inner_iterations',sum([s.history.inner_iterations]),'history',s.history,'preconditioner_graph_launches',s.preconditioner_graph_launches,'combined_array_bound_bytes',s.combined_array_bound_bytes,'maximum_solver_owned_device_bytes',s.maximum_solver_owned_device_bytes);
    end
end
