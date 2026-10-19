function r=checkMxOutputFailure(inputPath)
% Explicit allocation-site injection, not an opaque MATLAB OOM test.
inputs=load(inputPath);h=[];guard=onCleanup(@cleanup);
radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
h=radia_cuda_hmatrix('create',inputs.cases(1).config);
r=struct('accepted',false,'commands',struct('command',{},'allocation_sites',{},'failure_events',{}));
for name={'stats','download'}
    command=name{1};
    if strcmp(command,'download'),radia_cuda_hmatrix('upload',h,inputs.cases(1).initial);end
    radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
    value=radia_cuda_hmatrix(command,h); %#ok<NASGU>
    counts=radia_cuda_hmatrix('diagnostic_mx_allocation_report');
    count=counts.attempts;assert(count>0);events=0;
    for target=1:count
        for repeat=1:10
            radia_cuda_hmatrix('diagnostic_mx_allocation_failure',target);
            failed=false;
            try,value=radia_cuda_hmatrix(command,h);catch e
                assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
            end
            report=radia_cuda_hmatrix('diagnostic_mx_allocation_report');
            assert(failed && report.attempts==target);
            radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
            assert(radia_cuda_hmatrix('live_count')==1);events=events+1;
        end
    end
    value=radia_cuda_hmatrix(command,h); %#ok<NASGU>
    r.commands(end+1)=struct('command',command,'allocation_sites',count,'failure_events',events); %#ok<AGROW>
end
radia_cuda_hmatrix('destroy',h);h=[];
assert(radia_cuda_hmatrix('live_count')==0);r.accepted=true;
    function cleanup()
        radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
        if ~isempty(h),radia_cuda_hmatrix('destroy',h);end
    end
end
