function result=compareHDivAssembly(meshPath,referencePath,beforePath,afterPath,resultPath)
% Idle-100 diagnostic: actual FE setup, transfers, material solve and release.
inputs=load(referencePath);test=inputs.solver_cases(end);
result=struct('accepted',false,'cold_payloads',struct([]),'samples',struct([]));
try
    for k=1:2
        row=sample(k);if k==1,result.cold_payloads=row;else,result.cold_payloads(k)=row;end
    end
    order=[1 2 2 1 1 2 2 1];
    for k=1:numel(order)
        row=sample(order(k));if k==1,result.samples=row;else,result.samples(k)=row;end
    end
    assert(radia_cuda_hmatrix('live_count')==0);
    result.accepted=true;
catch e
    result.error=struct('identifier',e.identifier,'message',e.message);saveResult();rethrow(e)
end
saveResult();
    function r=sample(which)
        if exist('radia_cuda_hmatrix','file'),assert(radia_cuda_hmatrix('live_count')==0);end
        clear radia_cuda_hmatrix
        if contains(path,beforePath),rmpath(beforePath);end
        if contains(path,afterPath),rmpath(afterPath);end
        paths={beforePath,afterPath};addpath(paths{which});
        d=gpuDevice;freeBefore=d.AvailableMemory;
        solution=radia.gpu.solveHDivRT0(meshPath,[1 .2 .1],test.inv_chi,1e-11,100,'hacapk');
        freeAfter=d.AvailableMemory;
        coefficientError=norm(solution.coefficients-test.expected)/norm(test.expected);
        residual=norm(test.rhs-test.operator*solution.coefficients)/norm(test.rhs);
        assert(max(coefficientError,residual)<1e-10&&solution.true_relative_residual<=1e-11);
        assert(radia_cuda_hmatrix('live_count')==0);
        r=struct('mode',which,'whole_seconds',solution.whole_seconds,'coefficient_relative_error',coefficientError, ...
            'independent_true_residual',residual,'iterations',solution.iterations, ...
            'owned_peak_device_bytes',solution.combined_peak_owned_device_bytes, ...
            'available_device_bytes_before',freeBefore,'available_device_bytes_after',freeAfter);
    end
    function saveResult()
        file=fopen(resultPath,'w');assert(file>=0);closer=onCleanup(@()fclose(file)); %#ok<NASGU>
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
