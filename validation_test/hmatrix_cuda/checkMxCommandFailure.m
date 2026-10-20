function r=checkMxCommandFailure(inputPath,family,target)
% One enumerable output-fault attempt, for bounded/chunked external sweeps.
% target=0 enumerates this command only; positive targets must be reachable.
assert(isscalar(target) && target>=0 && fix(target)==target);
inputs=load(inputPath);h=[];guard=onCleanup(@cleanup);
assert(radia_cuda_hmatrix('live_count')==0);
for diagnostic={'c_allocation','host_allocation','device_allocation','resource','publication','mx_allocation'}
    radia_cuda_hmatrix(['diagnostic_' diagnostic{1} '_failure'],0);
end
probe=[];if isfield(inputs,'entry_probe'),probe=inputs.entry_probe;end
twoOutputs=false;info=struct();
switch family
    case 'live_count'
        args={family};
    case 'analytic_entry_batch'
        args={family,probe.cells,probe.faces,probe.rows,probe.columns};
    case 'analytic_aca_vector'
        args={family,probe.cells,probe.faces,[4;2;1;3],1,[0;1;0;0],0};
    case {'analytic_aca_compress','analytic_aca_leaf_compressed'}
        command='analytic_aca_compress';
        if strcmp(family,'analytic_aca_leaf_compressed'),command='analytic_aca_leaf';end
        args={command,probe.cells,probe.faces,2,1e-10,2};
    case 'analytic_aca_leaf_promoted'
        args={'analytic_aca_leaf',probe.cells,probe.faces,2,1e-10,1};
    case 'solve_linear_material'
        test=inputs.solver_cases(1);
        h=radia_cuda_hmatrix('create',test.config);
        args={family,h,test.mass,test.rhs,test.inv_chi,1e-11,100};
    case {'curved_entry_batch','curved_entry_batch_cached','curved_entry_batch_cached_info'}
        probe=inputs.curved_output_probe;
        assert(probe.cpu_operator_relative_error<=1e-10 && max(probe.setup_relative_errors)<=1e-11);
        command=family;
        if strcmp(family,'curved_entry_batch_cached_info')
            command='curved_entry_batch_cached';twoOutputs=true;
        end
        args={command,probe.config,probe.rows,probe.columns};
        if ~strcmp(command,'curved_entry_batch'),args{end+1}=probe.cache_budget_bytes;end
    otherwise
        error('Unsupported output allocation family: %s',family);
end
radia_cuda_hmatrix('diagnostic_mx_allocation_failure',target);
failed=false;identifier='';message='';value=[];
try
    [value,info]=invoke();
catch e
    failed=true;identifier=e.identifier;message=e.message;
end
% Capture before any recovery command can reset the command-local attempts.
report=radia_cuda_hmatrix('diagnostic_mx_allocation_report');
radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
if target==0
    assert(~failed && report.attempts>0);
else
    assert(failed && strcmp(identifier,'radia:cuda:hmatrix') && report.attempts==target, ...
        'Output injection target was not reached in the requested command');
    [value,info]=invoke();
end
numerical=validate(value);
if ~isempty(h),radia_cuda_hmatrix('destroy',h);h=[];end
r=struct('accepted',false,'family',family,'target',target,'failed',failed, ...
    'identifier',identifier,'message',message,'allocation_attempts',report.attempts, ...
    'numerical',numerical,'cache_info',info,'c',radia_cuda_hmatrix('diagnostic_c_allocation_report'), ...
    'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_report'), ...
    'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_report'), ...
    'resources',radia_cuda_hmatrix('diagnostic_resource_report'), ...
    'publication',radia_cuda_hmatrix('diagnostic_publication_report'), ...
    'live_count',radia_cuda_hmatrix('live_count'));
assert(r.live_count==0 && r.c.live_bytes==0 && r.cpp.live_bytes==0 && ...
    r.gpu.live_bytes==0 && r.gpu.free_errors==0 && ...
    r.resources.live_resources==0 && r.resources.free_errors==0 && ...
    r.publication.owned_output_arrays==0);
r.accepted=true;
    function n=validate(actual)
        n=struct('relative_error',0);
        switch family
            case 'live_count'
                expected=0;
            case {'curved_entry_batch','curved_entry_batch_cached','curved_entry_batch_cached_info'}
                expected=probe.expected;
                n.independent_mapped_operator_relative_error=probe.cpu_operator_relative_error;
                n.comparison_scope='Indexed/cached GPU entries; independent CPU FE operator checked in fixture';
                if twoOutputs
                    assert(info.cache_budget_bytes==probe.cache_budget_bytes && info.cache_tile_loads>=1);
                end
            case 'analytic_entry_batch'
                expected=probe.expected;
            case 'analytic_aca_vector'
                permutation=[4;2;1;3];expected=probe.aca_matrix(permutation(1),permutation).';
                expected(2)=0;
            case {'analytic_aca_compress','analytic_aca_leaf_compressed'}
                actual=actual.a2*actual.a1.';expected=probe.aca_matrix(1:2,3:4);
            case 'analytic_aca_leaf_promoted'
                assert(actual.promoted && actual.rank==0);
                actual=actual.dense;expected=probe.aca_matrix(1:2,3:4);
            case 'solve_linear_material'
                n.independent_true_residual=norm(test.rhs-test.operator*actual.x)/norm(test.rhs);
                n.solver_true_residual=actual.true_relative_residual;
                assert(n.independent_true_residual<=1e-10 && n.solver_true_residual<=1e-11);
                actual=actual.x;expected=test.expected;
        end
        n.relative_error=norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
        limit=1e-12;
        if contains(family,'compress'),limit=1e-8;end
        if strcmp(family,'solve_linear_material'),limit=1e-10;end
        assert(n.relative_error<=limit,'Declared numerical comparison rejected');
    end
    function [actual,cacheInfo]=invoke()
        cacheInfo=struct();
        if twoOutputs
            [actual,cacheInfo]=radia_cuda_hmatrix(args{:});
        else
            actual=radia_cuda_hmatrix(args{:});
        end
    end
    function cleanup()
        radia_cuda_hmatrix('diagnostic_mx_allocation_failure',0);
        if ~isempty(h),radia_cuda_hmatrix('destroy',h);end
    end
end
