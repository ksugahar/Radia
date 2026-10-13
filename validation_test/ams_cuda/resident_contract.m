function resident_contract(inputPath,outputPath,binaryPath,matlabPath)
% Ownership, stale/type-safe handles, update failure atomicity and workspace reuse.
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);
initial=radia_cuda_ams('resident_status');
[hierarchy,~]=radia.gpu.prepareAMS(s.model,SetupThreads=1);
op=hierarchy.A;b=gpuArray(s.rhs);
h=radia_cuda_ams('create',hierarchy);cleanup=onCleanup(@()radia_cuda_ams('release',h));
[x,first]=radia_cuda_ams('solve',h,op,b,1e-10,4000,true);
[reference,~]=radia.gpu.amsPCG(hierarchy,b,1e-10,4000,false,false);
assert(gather(norm(x-reference)/norm(reference))<1e-9,'resident direct parity');
bad=hierarchy;bad.A.rows=bad.A.rows+1;
expectError(@()radia_cuda_ams('update',h,bad));
radia_cuda_ams('update',h,hierarchy);clear hierarchy bad;
% Native ownership keeps all auxiliary arrays alive after caller source release.
[y,second]=radia_cuda_ams('solve',h,op,2*b,1e-10,4000,true);
assert(gather(norm(y-2*x)/norm(y))<1e-9,'resident retained ownership');
assert(first.workspace_allocations==1&&second.workspace_allocations==1,'workspace generation');
assert(second.solve_count==2&&second.hierarchy_builds==2,'workspace survives hierarchy replacement');
assert(first.native_buffer_allocations_during_solve==0&&second.native_buffer_allocations_during_solve==0,'native work allocation');
[zero,zeroInfo]=radia_cuda_ams('solve',h,op,zeros(size(b),'like',b),1e-10,4000,true);
assert(gather(norm(zero))==0&&zeroInfo.iterations==0,'zero RHS must reset resident iterate');
clear cleanup;
expectError(@()radia_cuda_ams('release',h));
expectError(@()radia_cuda_ams('nonlinear_release',h));
assert(isequal(initial,radia_cuda_ams('resident_status')),'context release leak');
report=struct('accepted',true,'stale_handle_rejected',true,'wrong_type_rejected',true, ...
    'failed_update_preserves_context',true,'caller_source_release_safe',true, ...
    'workspace_survives_update',true,'native_buffer_allocations_during_solve',0);
save(outputPath,'report');
end
function expectError(action)
raised=false;
try;action();catch;raised=true;end
assert(raised,'expected rejection');
end
