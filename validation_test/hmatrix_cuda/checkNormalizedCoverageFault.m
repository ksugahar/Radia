function result=checkNormalizedCoverageFault(root)
% Diagnostic build only. Enumerate the new site, then its frozen 5+10 class.
assert(~isfile(fullfile(root,'coverage-fault-samples.jsonl')));
[normal,fixture]=test_hmatrix_normalized_leaves();
assert(normal.accepted && normal.live_count==0);
device=gpuDevice;process=System.Diagnostics.Process.GetCurrentProcess();
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_resource_failure',0);
radia_cuda_hmatrix('diagnostic_publication_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0,'normalized.coverage');
h=radia_cuda_hmatrix('create_normalized',fixture.config);
counts=radia_cuda_hmatrix('diagnostic_host_allocation_report');
radia_cuda_hmatrix('destroy',h);
assert(counts.attempts==1 && strcmp(counts.site,'normalized.coverage'));
for target=1:counts.attempts,inject(target);end
% Prove the first unexecuted target is unreachable, rather than treating it as coverage.
radia_cuda_hmatrix('diagnostic_host_allocation_failure',counts.attempts+1,'normalized.coverage');
h=radia_cuda_hmatrix('create_normalized',fixture.config);
boundary=radia_cuda_hmatrix('diagnostic_host_allocation_report');
radia_cuda_hmatrix('destroy',h);assert(boundary.attempts==counts.attempts);
samples=zeros(15,2);
for iteration=1:15
    receipt=inject(1);
    wait(device);process.Refresh();
    samples(iteration,:)=[double(process.PrivateMemorySize64),double(device.TotalMemory-device.AvailableMemory)];
    record=struct('iteration',iteration,'warmup',iteration<=5,'sample',samples(iteration,:),'counters',receipt);
    fid=fopen(fullfile(root,'coverage-fault-samples.jsonl'),'a');assert(fid>=0);
    fprintf(fid,'%s\n',jsonencode(record));fclose(fid); % persist BEFORE judging
end
measured=samples(6:15,:);slopes=zeros(1,2);
for column=1:2,p=polyfit((1:10)',measured(:,column),1);slopes(column)=p(1);end
deltas=measured(end,:)-measured(1,:);
radia.gpu.writeHDivH5(fullfile(root,'coverage-fault-memory.h5'),struct('samples',samples), ...
    struct('units','SI','labels',{{}},'warmup',5,'measured',10,'slopes_bytes_per_run',slopes,'deltas_bytes',deltas));
assert(all(slopes<=2^20) && all(deltas<=2^20));
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
assert(radia_cuda_hmatrix('live_count')==0);
clear radia_cuda_hmatrix;
result=struct('accepted',true,'site','normalized.coverage','enumerated_fault_points',counts.attempts, ...
    'memory_fault_repeats',15,'warmup',5,'measured',10,'slopes_bytes_per_run',slopes, ...
    'deltas_bytes',deltas,'scope','New bitmap allocation only; not all-site faults or common100');

    function receipt=inject(target)
        radia_cuda_hmatrix('diagnostic_host_allocation_failure',target,'normalized.coverage');
        failed=false;identifier='';message='';
        try
            unexpected=radia_cuda_hmatrix('create_normalized',fixture.config);
            radia_cuda_hmatrix('destroy',unexpected);
        catch problem
            failed=true;identifier=problem.identifier;message=problem.message;
        end
        receipt=struct('host',radia_cuda_hmatrix('diagnostic_host_allocation_report'), ...
            'device',radia_cuda_hmatrix('diagnostic_device_allocation_report'), ...
            'resource',radia_cuda_hmatrix('diagnostic_resource_report'), ...
            'publication',radia_cuda_hmatrix('diagnostic_publication_report'));
        assert(failed && strcmp(identifier,'radia:cuda:hmatrix') && contains(lower(message),'alloc'));
        assert(receipt.host.attempts==target && strcmp(receipt.host.site,'normalized.coverage'));
        assert(receipt.host.live_bytes==0 && receipt.device.live_bytes==0 && receipt.device.free_errors==0);
        assert(receipt.resource.live_resources==0 && receipt.resource.free_errors==0);
        assert(receipt.publication.owned_output_arrays==0 && receipt.publication.live_handles==0);
        assert(radia_cuda_hmatrix('live_count')==0);
        radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
        operator=radia.gpu.HMatrixOperator(fixture.config);
        operator.upload(fixture.x);operator.apply();
        errorValue=norm(operator.download()-fixture.expected)/norm(fixture.expected);
        delete(operator);assert(errorValue<=1e-14);
        after=radia_cuda_hmatrix('diagnostic_host_allocation_report');assert(after.live_bytes==0);
    end
end
