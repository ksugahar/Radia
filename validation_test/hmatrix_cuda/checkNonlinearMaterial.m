function result=checkNonlinearMaterial(meshPath,cpuPath)
[cpu,~]=radia.gpu.readHDivH5(cpuPath);q=radia_cuda_projection('hdiv_material_quadrature',meshPath,2);
nq=numel(q.weights);values=reshape(q.evaluation*cpu.m(:),nq,3);weights=cpu.weights(:);
errors=struct('weights',norm(q.weights-weights)/norm(weights),'points',norm(q.points-cpu.points,'fro')/norm(cpu.points,'fro'),'values',norm(values-cpu.values,'fro')/norm(cpu.values,'fro'));
assert(max([errors.weights errors.points errors.values])<=1e-11);
law=radia.gpu.makeHDivBHMaterial(cpu.table);E=cell(1,3);for c=1:3,E{c}=gpuArray(q.evaluation((c-1)*nq+(1:nq),:));end
w=gpuArray(q.weights);out=radia.gpu.evaluateHDivBHMaterial(E,w,law,cpu.m(:),true);
errors.energy=abs(out.energy-cpu.energy)/max(abs(cpu.energy),realmin);errors.gradient=norm(out.gradient-cpu.gradient(:))/norm(cpu.gradient(:));errors.tangent=norm(out.tangent-cpu.tangent,'fro')/norm(cpu.tangent,'fro');
assert(max([errors.energy errors.gradient errors.tangent])<=1e-8);
law.device=struct('knots',gpuArray(law.knots),'coefficients',gpuArray(law.coefficients),'integral_base',gpuArray(law.integral_base));
cached=radia.gpu.evaluateHDivBHMaterial(E,w,law,cpu.m(:),true);
errors.device_cache=max([abs(cached.energy-out.energy)/max(abs(out.energy),realmin),norm(cached.gradient-out.gradient)/norm(out.gradient),norm(cached.tangent-out.tangent,'fro')/norm(out.tangent,'fro')]);
assert(errors.device_cache<=1e-13);

d=cpu.direction(:);h=1e-6*max(norm(cpu.m),1);plus=radia.gpu.evaluateHDivBHMaterial(E,w,law,cpu.m(:)+h*d,false);minus=radia.gpu.evaluateHDivBHMaterial(E,w,law,cpu.m(:)-h*d,false);
errors.energy_gradient=abs((plus.energy-minus.energy)/(2*h)-out.gradient'*d)/max(abs(out.gradient'*d),realmin);
errors.gradient_hessian=norm((plus.gradient-minus.gradient)/(2*h)-out.tangent*d)/max(norm(out.tangent*d),realmin);
assert(errors.energy_gradient<1e-5 && errors.gradient_hessian<1e-5);
result=struct('accepted',true,'dofs',q.dofs,'points',nq,'integration_order',q.integration_order,'errors',errors,'scope','Independent CPU material subgate only, no complete nonlinear solve acceptance');
end
