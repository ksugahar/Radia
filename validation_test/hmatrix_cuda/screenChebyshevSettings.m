function result=screenChebyshevSettings(mesh,cpu)
[ref,~]=radia.gpu.readHDivH5(cpu);settings=[8 .1;16 .1;16 .01;32 .01];rows={};
for j=1:size(settings,1)
    opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true,'FinalResidualOnly',true,'Preconditioner','mass-chebyshev','ChebyshevDegree',settings(j,1),'ChebyshevLowerFraction',settings(j,2)};
    whole=radia.gpu.solveNonlinearHDivElements(mesh,2,ref.applied_field(:)',ref.table,opts{:});verify(whole);
    obj=radia.gpu.PreparedHDivNonlinear(mesh,2,ref.applied_field(:)',ref.table,opts{:});owner=onCleanup(@()delete(obj));info=obj.stats();prepared=obj.solve();verify(prepared);
    warm=obj.solve(1,prepared.coefficients);verify(warm);assert(warm.nonlinear_iterations==0);
    clear owner obj;wait(gpuDevice);assert(radia_cuda_hmatrix('live_count')==0 && radia_cuda_projection('live_count')==0);
    rows{end+1}=struct('degree',settings(j,1),'lower_fraction',settings(j,2),'whole',whole,'prepared',prepared,'prepare',info,'initial_solution_reuse',warm);
end
result=struct('accepted',true,'cases',[rows{:}],'live_count',0,'scope','One whole and one zero-start prepared solve per configuration, screening only; no comparative timing acceptance');
    function verify(s)
        assert(norm(s.coefficients-ref.coefficients(:))/norm(ref.coefficients(:))<=1e-8 && s.nonlinear_relative_residual<=1e-12);
    end
end
