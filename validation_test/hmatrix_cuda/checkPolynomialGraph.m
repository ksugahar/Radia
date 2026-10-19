function result=checkPolynomialGraph(mesh,cpu,large)
% Same fixed SPD arithmetic, explicit graph execution and independent CPU gates.
[reference,~]=radia.gpu.readHDivH5(cpu);
field=reference.applied_field(:)';
if large,chunk=128;cache=ceil(1.6*2^30);tol=1e-12;else,chunk=7;cache=32*2^20;tol=1e-10;end
opts={'NonlinearTolerance',tol,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',cache,'MaterialChunkElements',chunk,'InexactNewton',true};
baseline=radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-polynomial');
graph=radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-polynomial-graph');
error=norm(graph.coefficients-baseline.coefficients)/norm(baseline.coefficients);
cpuError=norm(graph.coefficients-reference.coefficients(:))/norm(reference.coefficients(:));
assert(error<=1e-11 && cpuError<=1e-8 && graph.nonlinear_relative_residual<=tol);
assert(baseline.preconditioner_graph_launches==0 && graph.preconditioner_graph_launches>0);
obj=radia.gpu.PreparedHDivNonlinear(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-polynomial-graph');owner=onCleanup(@()delete(obj));
repeated=obj.solve();assert(norm(repeated.coefficients-graph.coefficients)/norm(graph.coefficients)<=1e-11);
rejected=false;try,obj.solve(1,ones(numel(reference.coefficients)+1,1));catch,rejected=true;end;assert(rejected);
recovered=obj.solve(1,repeated.coefficients);assert(recovered.nonlinear_iterations==0);
clear owner obj;wait(gpuDevice);assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'baseline',baseline,'graph',graph,'repeated',repeated,'graph_to_original_relative_error',error,'cpu_coefficient_relative_error',cpuError,'invalid_input_rejected',rejected,'live_count',radia_cuda_hmatrix('live_count'),'scope','Fixed eight-step SPD polynomial graph only, not full PCG graph; one-shot timings diagnostic','memory_limitations','Native buffer accounting excludes CUDA graph driver allocations and runtime/global peaks');
end
