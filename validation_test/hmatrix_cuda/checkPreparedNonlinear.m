function result=checkPreparedNonlinear(mesh,cpu,large)
% Independent CPU fixture comparison, explicit preparation/reuse/lifecycle gates.
[reference,~]=radia.gpu.readHDivH5(cpu);
if large
    field=[0 0 1e5];chunk=128;cache=ceil(1.6*2^30);preconditioner='mass-polynomial';tol=1e-12;
else
    field=[0 0 1e6];chunk=7;cache=32*2^20;preconditioner='mass-cholesky';tol=1e-10;
end
opts={'NonlinearTolerance',tol,'LinearTolerance',1e-12,'LinearMaxIterations',4000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',cache,'MaterialChunkElements',chunk,'Preconditioner',preconditioner};
device=gpuDevice;before=device.AvailableMemory;
one=radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:});
assert(one.live_count==0);validate(one);
cases=struct();
for inexact=[false true]
    timer=tic;obj=radia.gpu.PreparedHDivNonlinear(mesh,2,field,reference.table,opts{:},'InexactNewton',inexact);owner=onCleanup(@()delete(obj));
    info=obj.stats();first=obj.solve();validate(first);
    second=obj.solve();validate(second);
    assert(norm(first.coefficients-second.coefficients)/norm(first.coefficients)<=1e-12);
    retained=device.AvailableMemory;
    rejected=false;try,obj.solve(1,ones(numel(reference.coefficients)+1,1));catch,rejected=true;end;assert(rejected);
    % Invalid call must not corrupt immutable maps or native operator.
    afterFailure=obj.solve(1,second.coefficients);validate(afterFailure);
    assert(afterFailure.nonlinear_iterations==0);
    elapsed=toc(timer);clear owner obj;wait(device);assert(radia_cuda_hmatrix('live_count')==0);
    name='strict';if inexact,name='inexact';end
    cases.(name)=struct('prepare',info,'first',first,'second',second,'explicit_initial',afterFailure,'total_through_release_seconds',elapsed,'invalid_call_rejected',rejected,'available_device_memory_retained_sample',retained);
end
result=struct('accepted',true,'one_shot',one,'cases',cases,'available_device_memory_before_sample',before,'available_device_memory_after_release_sample',device.AvailableMemory,'live_count',radia_cuda_hmatrix('live_count'),'memory_limitations','Available-memory samples and controlled-array bounds, not global peak','scope','Same mesh and BH context repeated from zero; explicit initial solution tested separately; no process-cold claim');
    function validate(s)
        assert(norm(s.coefficients-reference.coefficients(:))/norm(reference.coefficients(:))<=1e-8,'Independent CPU coefficient mismatch');
        assert(norm(s.mean_magnetization-reference.mean(:))/norm(reference.mean(:))<=1e-8,'Independent CPU mean mismatch');
        assert(s.nonlinear_relative_residual<=tol);
        if ~isempty(s.history),assert(all([s.history.inner_true_residual]<=1e-6));end
    end
end
