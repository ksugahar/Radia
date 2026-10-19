function r=checkAllocationFailure(inputPath,kind,target,path,site)
% One bounded diagnostic case. Caller records process memory after each call.
if nargin<4,path='tiny-dense';end
if nargin<5,site='';end
input=load(inputPath);cfg=input.full_assembly_probe.config;
command='create_analytic_hacapk';test=[];numericError=[];
switch path
    case 'tiny-dense'
        cfg.cells=cfg.cells(1:24);cfg.faces=zeros(0,1);cfg.leaf_size=64;
    case 'adaptive'
        % Keep the independent fixture's separated geometry and existing ACA settings.
    case 'numeric-lowrank'
        index=find(strcmp({input.cases.name},'lowrank'),1);
        assert(~isempty(index),'Missing independent low-rank fixture');
        test=input.cases(index);cfg=test.config;command='create';
    otherwise,error('Unsupported allocation diagnostic path');
end
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_trapped_frame_exception',0);
radia_cuda_hmatrix('diagnostic_publication_failure',0);
if ~isempty(site),radia_cuda_hmatrix('diagnostic_host_allocation_failure',0,site);end
switch kind
    case 'c',radia_cuda_hmatrix('diagnostic_c_allocation_failure',target);
    case 'publication',radia_cuda_hmatrix('diagnostic_publication_failure',target);
    case 'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_failure',target);
    case 'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_failure',target,site);
    case 'matlab',radia_cuda_hmatrix('diagnostic_trapped_frame_exception',1);
    case 'baseline'
    otherwise,error('Unknown allocation failure class');
end
failed=false;identifier='';message='';h=[];
try,h=radia_cuda_hmatrix(command,cfg);
    if strcmp(kind,'baseline') && ~isempty(test)
        radia_cuda_hmatrix('upload',h,test.initial);
        radia_cuda_hmatrix('apply',h,1);
        actual=radia_cuda_hmatrix('download',h);
        numericError=norm(actual-test.expected_first)/norm(test.expected_first);
        assert(numericError<=1e-12,'Independent low-rank apply rejected');
    end
catch e,failed=true;identifier=e.identifier;message=e.message;
end
if ~isempty(h),radia_cuda_hmatrix('destroy',h);end
r=struct('path',path,'kind',kind,'target',target,'failed',failed,'identifier',identifier, ...
    'message',message,'allocation_site',site,'numeric_relative_error',numericError, ...
    'c',radia_cuda_hmatrix('diagnostic_c_allocation_report'), ...
    'cpp',radia_cuda_hmatrix('diagnostic_host_allocation_report'), ...
    'gpu',radia_cuda_hmatrix('diagnostic_device_allocation_report'), ...
    'live_count',radia_cuda_hmatrix('live_count'), ...
    'publication',radia_cuda_hmatrix('diagnostic_publication_report'));
radia_cuda_hmatrix('diagnostic_c_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_host_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_device_allocation_failure',0);
radia_cuda_hmatrix('diagnostic_trapped_frame_exception',0);
radia_cuda_hmatrix('diagnostic_publication_failure',0);
assert(r.live_count==0 && r.c.live_bytes==0 && r.cpp.live_bytes==0 && r.gpu.live_bytes==0 && r.gpu.free_errors==0 && r.publication.owned_output_arrays==0 && r.publication.live_handles==0);
if strcmp(kind,'baseline'),assert(~failed);
else
    assert(failed && strcmp(identifier,'radia:cuda:hmatrix'));
    if strcmp(kind,'publication'),assert(r.publication.reached_stage==target);end
    if strcmp(kind,'c'),assert(r.c.allocations==target && r.c.failure_code==1);end
    if strcmp(kind,'gpu'),assert(r.gpu.attempts==target);end
    if strcmp(kind,'cpp'),assert(r.cpp.attempts==target);end
    if strcmp(kind,'matlab'),assert(contains(message,'Injected trapped MATLAB exception'));end
end
r.accepted=true;
end
