function nonlinear_resident(inputPath,outputPath,binaryPath,matlabPath)
% Independent NGSolve-generated FE maps are validation inputs, not a GPU adapter.
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);d=gpuDevice();
before=d.AvailableMemory;
inputBytes=8*(numel(s.Cx)+numel(s.Cy)+numel(s.Cz)+numel(s.mass));
assert(10*inputBytes<before-max(2*1024^3,.2*d.TotalMemory), ...
    'radia:nonlinear:memory','Dense prototype exceeds conservative memory budget');
wait(d);uploadTimer=tic;
C={gpuArray(s.Cx),gpuArray(s.Cy),gpuArray(s.Cz)};w=gpuArray(s.weights);
M=gpuArray(s.mass);rhs=gpuArray(s.rhs);free=gpuArray(logical(s.model.free));
index=gpuArray(s.csrIndex);wait(d);uploadSeconds=toc(uploadTimer);
ownedInputBytes=8*(numel(s.Cx)+numel(s.Cy)+numel(s.Cz)+numel(s.mass)+numel(s.rhs)+numel(s.weights)+numel(s.csrIndex));
routes=[1 2 3 3 2 1];reports=cell(1,numel(routes));solutions=cell(1,numel(routes));
for sample=1:numel(routes)
    route=routes(sample);
    wait(d);total=tic;x=zeros(size(rhs),'like',rhs);hierarchy=[];history=[];
    builds=0;setupSeconds=0;solveSeconds=0;firstIterations=0;rebuildNext=false;
    for outer=1:30
        B={C{1}*x,C{2}*x,C{3}*x};b2=B{1}.^2+B{2}.^2+B{3}.^2;
        nu=s.nu0+s.gamma*b2;r=M*x-rhs;J=M;
        CB=B{1}.*C{1}+B{2}.*C{2}+B{3}.*C{3};
        for k=1:3
            r=r+C{k}'*(w.*nu.*B{k});
            J=J+C{k}'*((w.*nu).*C{k});
        end
        J=J+2*s.gamma*(CB'*(w.*CB));r(~free)=0;
        J(~free,:)=0;J(:,~free)=0;
        relative=gather(norm(r)/norm(rhs));
        fprintf('route=%d outer=%d residual=%.3e\n',route,outer,relative);
        if relative<=1e-9;break;end
        % Positive nu0, gamma and mass imply a positive consistent tangent.
        op=s.model.A;op.rowptr=gpuArray(op.rowptr);op.columns=gpuArray(op.columns);
        op.values=J(index);
        rebuild=isempty(hierarchy)||route==1||(route==3&&rebuildNext);
        if rebuild
            wait(d);t=tic;model=s.model;model.A.values=gather(op.values);
            [hierarchy,~]=radia.gpu.prepareAMS(model,SetupThreads=1);
            wait(d);setupSeconds=setupSeconds+toc(t);builds=builds+1;firstIterations=0;
        end
        wait(d);t=tic;[delta,info]=radia.gpu.amsPCG(hierarchy,-r,1e-10,2000,false,false,op);
        wait(d);solveSeconds=solveSeconds+toc(t);
        if firstIterations==0;firstIterations=info.iterations;end
        % Explicit bounded backtracking uses the physical energy on device.
        energy=@(z) magneticEnergy(z,C,w,M,rhs,s.nu0,s.gamma);
        old=energy(x);slope=gather(r'*delta);alpha=1;
        while energy(x+alpha*delta)>old+1e-4*alpha*slope
            alpha=alpha/2;
            assert(alpha>=2^-20,'radia:nonlinear:lineSearch','Line search failed');
        end
        x=x+alpha*delta;
        history=[history;outer relative info.iterations info.true_relative_residual alpha]; %#ok<AGROW>
    end
    assert(relative<=1e-9,'radia:nonlinear:convergence','Nonlinear residual failed');
    wait(d);solutions{sample}=gather(x);elapsed=toc(total);
    reports{sample}=struct('route',route,'sample',sample,'total_seconds',elapsed,'upload_seconds',uploadSeconds,'input_owned_bytes',ownedInputBytes,'observed_device_bytes',double(before-d.AvailableMemory),'setup_seconds',setupSeconds, ...
        'solve_seconds',solveSeconds,'builds',builds,'rebuild_rule','route3: next iteration rebuild if PCG exceeds max(30,2*first solve after rebuild)','history',history,'final_residual',relative);
end
save(outputPath,'solutions','reports');
end
function value=magneticEnergy(x,C,w,M,rhs,nu0,gamma)
b2=(C{1}*x).^2+(C{2}*x).^2+(C{3}*x).^2;
value=gather(sum(w.*(0.5*nu0*b2+0.25*gamma*b2.^2))+0.5*x'*M*x-rhs'*x);
end
