function record=checkHostSetupMemory(root,h,bh)
% Bounded validation only: fixed host setup, full shared Newton reuse/lifetime.
device=gpuDevice;wait(device);process=System.Diagnostics.Process.GetCurrentProcess();
processOwner=onCleanup(@()process.Dispose());
clock=tic;record=struct('accepted',false,'scope','Explicit host setup hex2 full Newton common100', ...
    'warm',5,'measured',100,'opaque_runtime_retention_zero',false,'cohorts',{{}});
solutions=cell(1,3);samples=cell(1,3);modes={'reference-diagonal','reference-chebyshev-graph','reference-diagonal'};
stream=fullfile(root,'memory-samples.jsonl');
assert(~isfile(stream));f=fopen(stream,'w');assert(f>=0);fclose(f);
try
    for cohort=1:3
        mode=modes{cohort};createClose=cohort==3;obj=[];
        if ~createClose,obj=make(mode);end
        owner=onCleanup(@()delete(obj));
        values=zeros(size(h.mass,1),100);observed=zeros(100,2);
        % Warm the common observer and append-only checkpoint before action warmup.
        for k=1:5,observerSample=measure();checkpoint(cohort,-k,observerSample);end
        for k=1:5,action();measure();end
        baseline=measure();
        for k=1:100
            assert(toc(clock)<195,'Frozen reserve exhausted');
            values(:,k)=action();observed(k,:)=measure();checkpoint(cohort,k,observed(k,:));
        end
        clear owner obj
        assert(radia_cuda_hmatrix('live_count')==0);
        solutions{cohort}=values;samples{cohort}=observed;
        tail=observed(end-24:end,:);t=(-12:12)';slopes=sum(t.*(tail-tail(1,:)),1)/sum(t.^2);delta=tail(end,:)-tail(1,:);
        item=struct('preconditioner',mode,'create_solve_destroy',createClose,'baseline',baseline, ...
            'final25_slope_bytes_per_run',slopes,'final25_delta_bytes',delta, ...
            'full100_delta_bytes',observed(end,:)-observed(1,:));
        record.cohorts{cohort}=item;
        % Persist samples and every solution BEFORE judging a gate.
        persistSnapshot(sprintf('cohort-%d',cohort));
        assert(all(slopes<=65536) && all(delta<=1048576),'Frozen memory gate failed');
    end
    record.accepted=true;
catch ex
    record.error=getReport(ex,'extended');
    persistSnapshot('failure');
    rethrow(ex);
end
record.elapsed_seconds=toc(clock);persistSnapshot('final');
f=fopen(fullfile(root,'memory-result.json'),'w');assert(f>=0);fprintf(f,'%s',jsonencode(record));fclose(f);
    function persistSnapshot(tag)
        fields=struct('protocol_version',uint32(1));
        for c=1:numel(solutions)
            if ~isempty(solutions{c})
                fields.(sprintf('solutions_%d',c))=solutions{c};
                fields.(sprintf('samples_%d',c))=samples{c};
            end
        end
        % Immutable checkpoint, including failed-gate arrays; never overwrite.
        radia.gpu.writeHDivH5(string(fullfile(root,['memory-' tag '.h5'])),fields, ...
            struct('units','SI','labels',{{'body','exterior'}},'record',record));
    end
    function prepared=make(mode)
        prepared=radia.gpu.PreparedHDivNonlinear('',2,[1000 -2000 5000],bh, ...
            Assembly='ngsolve-host-normalized',HostSetup=h,Preconditioner=mode, ...
            NonlinearTolerance=5e-14,LinearTolerance=5e-14);
    end
    function x=action()
        if createClose
            temporary=make(mode);release=onCleanup(@()delete(temporary));v=temporary.solve();
        else
            v=obj.solve();
        end
        assert(v.nonlinear_relative_residual<=1e-11);
        x=v.coefficients;
    end
    function value=measure()
        wait(device);process.Refresh();value=[double(process.PrivateMemorySize64),double(device.TotalMemory-device.AvailableMemory)];
    end
    function checkpoint(c,i,value)
        f=fopen(stream,'a');assert(f>=0);closeFile=onCleanup(@()fclose(f));
        fprintf(f,'%s\n',jsonencode(struct('cohort',c,'iteration',i,'memory',value,'live_count',radia_cuda_hmatrix('live_count'))));
    end
end
