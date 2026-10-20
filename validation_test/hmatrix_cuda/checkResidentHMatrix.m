function result=checkResidentHMatrix(inputPath,resultPath,progressPath)
% Numerical/lifecycle subgate only, no assembly or performance claim.
arguments
    inputPath (1,:) char
    resultPath (1,:) char
    progressPath (1,:) char = ''
end
sequence=0;
inputs=load(inputPath);result=struct('accepted',false,'cases',struct([]));
try
    if isfield(inputs,'full_assembly_probe')
        test=inputs.full_assembly_probe;
        failed=false;try,unused=native('create_analytic_hacapk',[test.config test.config]);catch e %#ok<NASGU>
            assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
        end
        assert(failed&&native('live_count')==0,'Struct-array configuration must be rejected');
        invalidBudgets=0;
        for value=[0 -1 NaN 2*2^30+1 .5]
            invalid=test.config;invalid.frame_coefficient_budget_bytes=value;failed=false;
            try,unused=native('create_analytic_hacapk',invalid);catch e
                assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
            end
            assert(failed && native('live_count')==0);invalidBudgets=invalidBudgets+1;
        end
        lowBudget=test.config;lowBudget.frame_coefficient_budget_bytes=1;
        for attempt=1:3
            failed=false;
            try,unused=native('create_analytic_hacapk',lowBudget);catch e
                assert(strcmp(e.identifier,'radia:cuda:hmatrix') && contains(e.message,'HACApK coefficient frame allowance exceeded'));
                failed=true;
            end
            assert(failed && native('live_count')==0,'Rejected frame must not publish a handle');
        end
        result.frame_budget_rejection=struct('attempts',3,'invalid_budget_values_rejected',invalidBudgets,'live_count_after_rejection',0,'scope','Real module TLS callback and C frame unwind, not pre-parser rejection');
        for value=[0 -1 NaN 2*2^30+1 .5]
            invalid=test.config;invalid.frame_working_payload_budget_bytes=value;failed=false;
            try,unused=native('create_analytic_hacapk',invalid);catch e
                assert(strcmp(e.identifier,'radia:cuda:hmatrix'));failed=true;
            end
            assert(failed&&native('live_count')==0);
        end
        lowWorking=test.config;lowWorking.frame_working_payload_budget_bytes=1;
        for attempt=1:3
            failed=false;
            try,unused=native('create_analytic_hacapk',lowWorking);catch e
                assert(strcmp(e.identifier,'radia:cuda:hmatrix')&&contains(e.message,'HACApK coefficient/ACA working payload allowance exceeded'));
                failed=true;
            end
            assert(failed&&native('live_count')==0);
        end
        result.frame_working_payload_rejection=struct('attempts',3,'invalid_budget_values_rejected',5,'live_count_after_rejection',0,'scope','Actual C frame admission before fill; coefficient allowance remains default; subsequent normal create/apply must recover');
        h=native('create_analytic_hacapk',test.config);
        free=onCleanup(@()native('destroy',h));
        native('upload',h,test.initial);
        native('apply',h,1);actual=native('download',h);
        difference=norm(actual-test.expected)/max(norm(test.expected),realmin);
        assert(difference<=1e-8,'Full CUDA HACApK assembly differs from CPU direct entries');
        result.full_assembly=struct('relative_difference',difference,'stats',native('stats',h),'scope','Constant affine tetra charges; full native hierarchy and CUDA entries; not general HDiv physics');
        if isfield(test,'expected_matrix')
            n=size(test.expected_matrix,1);assembled=zeros(n,n);
            for column=1:n
                basis=zeros(n,1);basis(column)=1;
                native('upload',h,basis);
                native('apply',h,1);
                assembled(:,column)=native('download',h);
            end
            operatorError=norm(assembled-test.expected_matrix,'fro')/max(norm(test.expected_matrix,'fro'),realmin);
            result.full_assembly.operator_frobenius_relative_error=operatorError;
            result.full_assembly.operator_columns_checked=n;
            result.full_assembly.symmetry_relative_error=norm(assembled-assembled.','fro')/max(norm(assembled,'fro'),realmin);
            assert(operatorError<=1e-8,'Full CUDA operator differs from independent CPU direct matrix');
        end
        clear free
        assert(native('live_count')==0);
    end
    for k=1:numel(inputs.cases)
        test=inputs.cases(k);h=create(test.config);
        free=onCleanup(@()native('destroy',h));
        mustFail(@()native('apply',h,1));
        mustFail(@()native('upload',h,nan(size(test.initial))));
        native('upload',h,test.initial);
        native('apply',h,1);
        first=native('download',h);
        mustFail(@()native('upload',h,nan(size(test.initial))));
        unchanged=native('download',h);
        assert(isequal(first,unchanged));
        mustFail(@()native('apply',h,0));
        mustFail(@()native('apply',uint64(1),1));
        mustFail(@()native('apply',double(h),1));
        native('apply',h,2);
        third=native('download',h);
        firstError=norm(first-test.expected_first)/max(norm(test.expected_first),realmin);
        thirdError=norm(third-test.expected_third)/max(norm(test.expected_third),realmin);
        assert(max(firstError,thirdError)<=1e-12,'CPU/CUDA leaf application mismatch');
        stats=native('stats',h);
        launches=12;if strcmp(test.config.operator_kind,'demag'),launches=18;end
        assert(stats.applications==3&&stats.kernel_launches==launches&&stats.host_synchronizations_per_apply==0);
        caseResult=struct('name',test.name,'first_relative_error',firstError, ...
            'third_relative_error',thirdError,'stats',stats,'invalid_input_preserved_state',true);
        if k==1,result.cases=caseResult;else,result.cases(k)=caseResult;end
        clear free
        mustFail(@()native('apply',h,1));
        mustFail(@()native('destroy',h));
    end
    if isfield(inputs,'solver_cases')
        result.material_solver=struct([]);
        for k=1:numel(inputs.solver_cases)
            test=inputs.solver_cases(k);h=create(test.config);free=onCleanup(@()native('destroy',h));
            solution=native('solve_linear_material',h,test.mass,test.rhs,test.inv_chi,1e-11,100);
            error=norm(solution.x-test.expected)/norm(test.expected);
            residual=norm(test.rhs-test.operator*solution.x)/norm(test.rhs);
            assert(error<=1e-10&&residual<=1e-10&&solution.true_relative_residual<=1e-11);
            zero=native('solve_linear_material',h,test.mass,zeros(size(test.rhs)),test.inv_chi,1e-11,100);
            assert(norm(zero.x)==0&&zero.iterations==0);
            mustFail(@()native('solve_linear_material',h,test.mass,test.rhs,test.inv_chi,1e-11,1));
            native('upload',h,test.rhs);before=native('download',h);
            bad=test.mass;bad(1,2)=bad(1,2)+1;
            mustFail(@()native('solve_linear_material',h,bad,test.rhs,test.inv_chi,1e-11,100));
            assert(isequal(before,native('download',h)));
            row=struct('relative_error',error,'independent_true_residual',residual,'solver_true_residual',solution.true_relative_residual, ...
                'iterations',solution.iterations,'solver_owned_device_bytes',solution.solver_owned_device_bytes, ...
                'scope',test.scope,'mean_magnetization_relative_error',norm(test.mean_operator*(solution.x-test.expected))/max(norm(test.mean_operator*test.expected),realmin));
            if k==1,result.material_solver=row;else,result.material_solver(k)=row;end
            clear free
        end
    end
    probe=inputs.entry_probe;
    entries=native('analytic_entry_batch',probe.cells,probe.faces,probe.rows,probe.columns);
    entryError=norm(entries-probe.expected)/max(norm(probe.expected),realmin);
    assert(entryError<=1e-12,'Indexed analytic entry parity failed');
    result.indexed_entries=struct('count',numel(entries),'relative_error',entryError,'workspace_capacity',256);
    permutation=[4;2;1;3];mask=[0;1;0;0];acaError=0;
    for direction=0:1
        for pivot=1:4
            actual=native('analytic_aca_vector',probe.cells,probe.faces,permutation,pivot,mask,direction);
            if direction==0,expected=probe.aca_matrix(permutation(pivot),permutation).';
            else,expected=probe.aca_matrix(permutation,permutation(pivot));end
            expected(mask~=0)=0;
            acaError=max(acaError,norm(actual-expected)/max(norm(expected),realmin));
        end
    end
    assert(acaError<=1e-12,'ACA C dispatch to CUDA parity failed');
    result.aca_cuda_vector=struct('cases',8,'relative_error',acaError,'scope','Active vector C dispatch only, not compressed ACA matrix');
    factors=native('analytic_aca_compress',probe.cells,probe.faces,2,1e-10,2);
    reconstructed=factors.a2*factors.a1.';expected=probe.aca_matrix(1:2,3:4);
    compressionError=norm(reconstructed-expected,'fro')/max(norm(expected,'fro'),realmin);
    result.aca_cuda_compression=struct('rank',factors.rank,'relative_error',compressionError, ...
        'entry_batches',factors.entry_batches,'entry_count',factors.entry_count, ...
        'owned_device_bytes',factors.owned_device_bytes,'scope','One raw physical offdiagonal block, not full adaptive matrix or solver');
    assert(compressionError<=1e-8,'Existing ACA+ with CUDA entry compression parity failed');
    if isfield(inputs,'aca_large_probe')
        probeLarge=inputs.aca_large_probe;
        factors=native('analytic_aca_compress',probeLarge.cells,probeLarge.faces,probeLarge.split,1e-10,min(32,probeLarge.split));
        relativeError=norm(factors.a2*factors.a1.'-probeLarge.expected,'fro')/norm(probeLarge.expected,'fro');
        result.admissible_block=struct('rows',probeLarge.split,'columns',probeLarge.split,'rank',factors.rank, ...
            'relative_error',relativeError,'entry_count',factors.entry_count,'entry_batches',factors.entry_batches, ...
            'owned_device_bytes',factors.owned_device_bytes,'scope','Separated constant-tet block above default 65-charge leaf threshold', ...
            'cpu_leaf_stats',probeLarge.cpu_leaf_stats);
        assert(relativeError<=1e-8,'Admissible ACA block physical CPU parity failed');
        leaf=native('analytic_aca_leaf',probeLarge.cells,probeLarge.faces,probeLarge.split,1e-10,32);
        assert(leaf.promoted&&leaf.rank==0&&leaf.aca_rank==32);
        denseError=norm(leaf.dense-probeLarge.expected,'fro')/norm(probeLarge.expected,'fro');
        [rr,cc]=ndgrid(1:probeLarge.split,probeLarge.split+(1:probeLarge.split));
        direct=native('analytic_entry_batch',probeLarge.cells,probeLarge.faces,rr(:),cc(:));
        directError=norm(leaf.dense(:)-direct)/max(norm(direct),realmin);

        result.dense_promotion=struct('relative_error',denseError,'direct_cuda_relative_error',directError,'cpu_physical_limit',1e-8,'aca_rank',leaf.aca_rank, ...
            'entry_count',leaf.entry_count,'entry_batches',leaf.entry_batches, ...
            'owned_device_bytes',leaf.owned_device_bytes,'scope','Exact CUDA leaf; full adaptive assembly pending');
        assert(directError<=1e-12,'Dense promotion changed direct CUDA entries');
        assert(denseError<=1e-8,'Promoted leaf physical CPU parity failed');

    end
    assert(native('live_count')==0);
    stale=h;
    clear radia_cuda_hmatrix
    h=native('create',inputs.cases(1).config);
    free=onCleanup(@()native('destroy',h));
    assert(h~=stale);mustFail(@()native('apply',stale,1));
    clear free
    assert(native('live_count')==0);
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
    assert(native('live_count')==0);
    operator=radia.gpu.HMatrixOperator(inputs.cases(4).config);
    free=onCleanup(@()delete(operator));operator.upload(inputs.cases(4).initial);operator.apply(3);
    wrapper=operator.download();assert(norm(wrapper-inputs.cases(4).expected_third)<=1e-12*norm(inputs.cases(4).expected_third));
    clear free
    assert(native('live_count')==0);
    bad=inputs.cases(end).config;bad.cells(:)=0;bad.faces=zeros(0,1);mustFailCreate(bad);
    bad=inputs.cases(end).config;bad.image_masks=1;mustFailCreate(bad);
    assert(native('live_count')==0);
    result.wrapper_checked=true;result.accepted=true;result.live_count=0;result.reload_stale_handle_rejected=true;
    result.invalid_create_released=true;
catch failure
    result.error=struct('identifier',failure.identifier,'message',failure.message);
    saveResult();rethrow(failure)
end
saveResult();
    function h=create(config)
        command='create';if isfield(config,'assembly'),command='create_analytic_constant_block';end
        h=native(command,config);
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
    function varargout=native(varargin)
        sequence=sequence+1; current=sequence;
        mark(current,varargin{1},'enter','');
        try
            [varargout{1:nargout}]=radia_cuda_hmatrix(varargin{:});
        catch failure
            mark(current,varargin{1},'error',failure.identifier);
            rethrow(failure)
        end
        mark(current,varargin{1},'return','');
    end
    function mark(current,command,state,identifier)
        if isempty(progressPath),return;end
        % Close each append before entering native code: a killed worker leaves
        % its last command durable without requiring an Engine response.
        file=fopen(progressPath,'a');assert(file>=0,'Cannot append progress');
        closer=onCleanup(@()fclose(file)); %#ok<NASGU>
        fprintf(file,'%s\n',jsonencode(struct('sequence',current, ...
            'command',command,'state',state,'identifier',identifier)));
    end
    function saveResult()
        file=fopen(resultPath,'w');assert(file>=0);closer=onCleanup(@()fclose(file)); %#ok<NASGU>
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
