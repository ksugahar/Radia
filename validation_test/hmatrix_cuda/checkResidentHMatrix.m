function result=checkResidentHMatrix(inputPath,resultPath)
% Numerical/lifecycle subgate only, no assembly or performance claim.
arguments
    inputPath (1,:) char
    resultPath (1,:) char
end
inputs=load(inputPath);result=struct('accepted',false,'cases',struct([]));
try
    for k=1:numel(inputs.cases)
        test=inputs.cases(k);h=create(test.config);
        free=onCleanup(@()radia_cuda_hmatrix('destroy',h));
        mustFail(@()radia_cuda_hmatrix('apply',h,1));
        mustFail(@()radia_cuda_hmatrix('upload',h,nan(size(test.initial))));
        radia_cuda_hmatrix('upload',h,test.initial);
        radia_cuda_hmatrix('apply',h,1);
        first=radia_cuda_hmatrix('download',h);
        mustFail(@()radia_cuda_hmatrix('upload',h,nan(size(test.initial))));
        unchanged=radia_cuda_hmatrix('download',h);
        assert(isequal(first,unchanged));
        mustFail(@()radia_cuda_hmatrix('apply',h,0));
        mustFail(@()radia_cuda_hmatrix('apply',uint64(1),1));
        mustFail(@()radia_cuda_hmatrix('apply',double(h),1));
        radia_cuda_hmatrix('apply',h,2);
        third=radia_cuda_hmatrix('download',h);
        firstError=norm(first-test.expected_first)/max(norm(test.expected_first),realmin);
        thirdError=norm(third-test.expected_third)/max(norm(test.expected_third),realmin);
        assert(max(firstError,thirdError)<=1e-12,'CPU/CUDA leaf application mismatch');
        stats=radia_cuda_hmatrix('stats',h);
        launches=12;if strcmp(test.config.operator_kind,'demag'),launches=18;end
        assert(stats.applications==3&&stats.kernel_launches==launches&&stats.host_synchronizations_per_apply==0);
        caseResult=struct('name',test.name,'first_relative_error',firstError, ...
            'third_relative_error',thirdError,'stats',stats,'invalid_input_preserved_state',true);
        if k==1,result.cases=caseResult;else,result.cases(k)=caseResult;end
        clear free
        mustFail(@()radia_cuda_hmatrix('apply',h,1));
        mustFail(@()radia_cuda_hmatrix('destroy',h));
    end
    assert(radia_cuda_hmatrix('live_count')==0);
    stale=h;
    clear radia_cuda_hmatrix
    h=radia_cuda_hmatrix('create',inputs.cases(1).config);
    free=onCleanup(@()radia_cuda_hmatrix('destroy',h));
    assert(h~=stale);mustFail(@()radia_cuda_hmatrix('apply',stale,1));
    clear free
    assert(radia_cuda_hmatrix('live_count')==0);
    mustFailCreate(1);
    mustFailCreate(repmat(inputs.cases(1).config,1,2));
    bad=inputs.cases(1).config;bad.permutation(:)=1;
    mustFailCreate(bad);
    bad=inputs.cases(1).config;bad.scale(1)=0;
    mustFailCreate(bad);
    bad=inputs.cases(1).config;bad.leaves(1).row=bad.dimension+1;
    mustFailCreate(bad);
    bad=inputs.cases(1).config;bad.leaves(1).a1=nan;
    mustFailCreate(bad);
    assert(radia_cuda_hmatrix('live_count')==0);
    operator=radia.gpu.HMatrixOperator(inputs.cases(4).config);
    free=onCleanup(@()delete(operator));operator.upload(inputs.cases(4).initial);operator.apply(3);
    wrapper=operator.download();assert(norm(wrapper-inputs.cases(4).expected_third)<=1e-12*norm(inputs.cases(4).expected_third));
    clear free
    assert(radia_cuda_hmatrix('live_count')==0);
    bad=inputs.cases(end).config;bad.cells(:)=0;bad.faces=zeros(0,1);mustFailCreate(bad);
    bad=inputs.cases(end).config;bad.image_masks=1;mustFailCreate(bad);
    assert(radia_cuda_hmatrix('live_count')==0);
    result.wrapper_checked=true;result.accepted=true;result.live_count=0;result.reload_stale_handle_rejected=true;
    result.invalid_create_released=true;
catch failure
    result.error=struct('identifier',failure.identifier,'message',failure.message);
    saveResult();rethrow(failure)
end
saveResult();
    function h=create(config)
        command='create';if isfield(config,'assembly'),command='create_analytic_constant_block';end
        h=radia_cuda_hmatrix(command,config);
    end
    function mustFail(call)
        failed=false;try,call();catch e
            assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
        end
        assert(failed,'Expected a checked native rejection');
    end
    function mustFailCreate(config)
        failed=false;try,unused=create(config);catch e %#ok<NASGU>
            assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
        end
        assert(failed,'Expected a checked create rejection');
    end
    function saveResult()
        file=fopen(resultPath,'w');assert(file>=0);closer=onCleanup(@()fclose(file)); %#ok<NASGU>
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
