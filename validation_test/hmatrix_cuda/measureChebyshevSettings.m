function result=measureChebyshevSettings(mesh,cpu,degree,lower,compareGraph)
if nargin<5,compareGraph=false;end
assert(islogical(compareGraph)&&isscalar(compareGraph));
[ref,~]=radia.gpu.readHDivH5(cpu);settings=[8 .1;degree lower];device=gpuDevice;before=device.AvailableMemory;samples={};first={};
preconditioners={'mass-chebyshev','mass-chebyshev'};
if compareGraph,settings=[degree lower;degree lower];preconditioners{2}='mass-chebyshev-graph';end
for mode=1:2,first{mode}=runWhole(mode);end
for cycle=1:2
    for mode=[1 2 2 1]
        row=runWhole(mode);row.mode=mode;row.cycle=cycle;samples{end+1}=row;
    end
end
samples=[samples{:}];series={};
for mode=1:2
    timer=tic;opts=options(mode);obj=radia.gpu.PreparedHDivNonlinear(mesh,2,ref.applied_field(:)',ref.table,opts{:});owner=onCleanup(@()delete(obj));info=obj.stats();runs={};
    for k=1:3,runs{end+1}=record(obj.solve());end
    release=tic;clear owner obj;wait(device);releaseSeconds=toc(release);assert(radia_cuda_hmatrix('live_count')==0);
    series{mode}=struct('prepare',info,'runs',[runs{:}],'release_seconds',releaseSeconds,'construction_three_solves_release_seconds',toc(timer));
end
result=struct('accepted',true,'settings',settings,'compare_graph',compareGraph,'preconditioners',{preconditioners},'first_protocol_samples',[first{:}],'samples',samples,'baseline_whole_median_seconds',median([samples([samples.mode]==1).whole_seconds]),'candidate_whole_median_seconds',median([samples([samples.mode]==2).whole_seconds]),'prepared_series',[series{:}],'available_device_memory_before_sample',before,'available_device_memory_after_sample',device.AvailableMemory,'live_count',radia_cuda_hmatrix('live_count'),'scope','Same binary preconditioner comparison; first calls plus two warm ABBA cycles; whole FE/HMatrix/transfers/nonlinear solve/recovery/release; prepared three zero-start solves exclude construction/release; no CPU speedup or process cold','memory_limitations','Controlled native/MATLAB bounds and device samples, not host/global/library/graph-driver peak');
    function opts=options(mode)
        opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true,'FinalResidualOnly',true,'Preconditioner',preconditioners{mode},'ChebyshevDegree',settings(mode,1),'ChebyshevLowerFraction',settings(mode,2)};
    end
    function row=runWhole(mode)
        opts=options(mode);row=record(radia.gpu.solveNonlinearHDivElements(mesh,2,ref.applied_field(:)',ref.table,opts{:}));assert(radia_cuda_hmatrix('live_count')==0);
    end
    function row=record(s)
        e=norm(s.coefficients-ref.coefficients(:))/norm(ref.coefficients(:));assert(e<=1e-8 && s.nonlinear_relative_residual<=1e-12);
        if strcmp(s.preconditioner,'mass-chebyshev-graph'),assert(s.preconditioner_graph_launches>0);else,assert(s.preconditioner_graph_launches==0);end
        row=struct('coefficients',s.coefficients,'mean_magnetization',s.mean_magnetization,'preconditioner',s.preconditioner,'graph_launches',s.preconditioner_graph_launches,'degree',s.chebyshev_degree,'lower_fraction',s.chebyshev_lower_fraction,'whole_seconds',s.whole_seconds,'nonlinear_updates',s.nonlinear_iterations,'inner_linear_iterations',sum([s.history.inner_iterations]),'cpu_coefficient_relative_error',e,'native_nonlinear_relative_residual',s.nonlinear_relative_residual,'combined_array_bound_bytes',s.combined_array_bound_bytes,'solver_owned_device_bytes',s.maximum_solver_owned_device_bytes);
    end
end
