function r=checkResourceFailure(inputPath,preconditioner,target,kind,site)
% Complete create/solve/destroy unit; CPU operator/solution are independent fixtures.
% Optional kind/site enumerate host/device as well as graph resource failures.
% Configure before create: resetting counters with live operators is forbidden.
if nargin<4,kind='resource';end
if nargin<5,site='';end
assert(any(strcmp(kind,{'resource','cpp','gpu','c'})));
assert(isempty(site)||strcmp(kind,'cpp'));
input=load(inputPath);test=input.solver_cases(1);
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_resource_failure',0);
switch kind
    case 'resource',radia_cuda_hmatrix('diagnostic_resource_failure',target);
    case 'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_failure',target,site);
    case 'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_failure',target);
    case 'c',radia_cuda_hmatrix('diagnostic_c_allocation_failure',target);
end
radia_cuda_hmatrix('diagnostic_trapped_frame_exception',0);
h=[];failed=false;identifier='';message='';errorValue=[];residual=[];createCounts=struct();
try
    command='create';if isfield(test.config,'assembly'),command='create_analytic_constant_block';end
    h=radia_cuda_hmatrix(command,test.config);
    createCounts=struct('c',radia_cuda_hmatrix('diagnostic_c_allocation_report'), ...
        'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_report'), ...
        'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_report'), ...
        'resources',radia_cuda_hmatrix('diagnostic_resource_report'));
    solution=radia_cuda_hmatrix('solve_linear_material',h,test.mass,test.rhs,test.inv_chi,1e-11,200,preconditioner);
    errorValue=norm(solution.x-test.expected)/norm(test.expected);
    residual=norm(test.rhs-test.operator*solution.x)/norm(test.rhs);
    assert(errorValue<=1e-10 && residual<=1e-10 && solution.true_relative_residual<=1e-11);
catch e
    failed=true;identifier=e.identifier;message=e.message;
end
if ~isempty(h),radia_cuda_hmatrix('destroy',h);end
report=radia_cuda_hmatrix('diagnostic_resource_report');
gpu=radia_cuda_hmatrix('diagnostic_device_allocation_report');
cReport=radia_cuda_hmatrix('diagnostic_c_allocation_report');
cpp=radia_cuda_hmatrix('diagnostic_host_allocation_report');
live=radia_cuda_hmatrix('live_count');
assert(live==0 && report.live_resources==0 && report.free_errors==0 && gpu.live_bytes==0 && gpu.free_errors==0 && cReport.live_bytes==0 && cpp.live_bytes==0);
switch kind
    case 'resource',attempts=report.attempts;
    case 'cpp',attempts=cpp.attempts;assert(strcmp(cpp.filter,site));
    case 'gpu',attempts=gpu.attempts;
    case 'c',attempts=cReport.allocations;
end
if target==0,assert(~failed);else,assert(failed && strcmp(identifier,'radia:cuda:hmatrix') && attempts==target);end
r=struct('accepted',true,'preconditioner',preconditioner,'target',target,'failed',failed, ...
    'identifier',identifier,'message',message,'resources',report,'gpu',gpu,'live_count',live, ...
    'cpu_error',errorValue,'independent_true_residual',residual, ...
    'kind',kind,'allocation_site',site,'allocation_attempts',attempts, ...
    'create_counts',createCounts,'c',cReport,'cpp',cpp);
radia_cuda_hmatrix('diagnostic_resource_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
end
