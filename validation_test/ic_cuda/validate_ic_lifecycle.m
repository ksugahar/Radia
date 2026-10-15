function result=validate_ic_lifecycle(simulateMismatch)
% Test-only build simulates selected-device rejection on single-GPU hosts.
if nargin==0;simulateMismatch=false;end
assert(radia_cuda_ic('live_count')==0,'Lifecycle test requires an owned empty runtime');
dev=gpuDevice(); original=dev.Index;
op=struct('rows',2.,'cols',2.,'rowptr',gpuArray(int32([0;1;2])), ...
    'columns',gpuArray(int32([0;1])),'values',gpuArray([2.;3.]));
[h,~]=radia_cuda_ic('create',op,1.05);
path='same-device'; caller=original;
if simulateMismatch
    radia_cuda_ic('test_device_mismatch',h);path='test-only selected-device mismatch simulation';
elseif gpuDeviceCount('available')>1
    available=gpuDeviceTable(); caller=available.Index(find(available.Index~=original & available.DeviceAvailable,1));
    gpuDevice(caller);path='physical device switch';
end
restore=onCleanup(@()gpuDevice(original));
if simulateMismatch || caller~=original
    rhs=gpuArray([1.;1.]);
    mustReject(@()apply(h,rhs));mustReject(@()solve(h,rhs));
end
radia_cuda_ic('release',h);
assert(gpuDevice().Index==caller,'Release changed the caller device');
assert(radia_cuda_ic('live_count')==0);
mustReject(@()radia_cuda_ic('release',h));
mustReject(@()radia_cuda_ic('release',uint64(1)));
mustReject(@()radia_cuda_ic('release',double(h)));
% Lock balance is observable: clear must unload the native module.
clear radia_cuda_ic;
[~,loaded]=inmem('-completenames');
assert(~any(contains(loaded,'radia_cuda_ic.mex')),'IC MEX remains locked');
clear restore;
% A stale generation must not revive after unloading/reloading the MEX.
op.values=gpuArray([2.;3.]);op.rowptr=gpuArray(int32([0;1;2]));op.columns=gpuArray(int32([0;1]));
[fresh,~]=radia_cuda_ic('create',op,1.05);
assert(fresh~=h);mustReject(@()radia_cuda_ic('release',h));radia_cuda_ic('release',fresh);
assert(radia_cuda_ic('live_count')==0);clear radia_cuda_ic;
result=struct('passed',true,'path',path,'caller_device_restored',true, ...
    'live_count',0,'mex_unloaded',true,'stale_rejected',true);
end
function apply(h,rhs)
x=radia_cuda_ic('apply',h,rhs); %#ok<NASGU>
end
function solve(h,rhs)
[x,info]=radia_cuda_ic('solve',h,rhs,1e-10,100); %#ok<ASGLU>
end
function mustReject(f)
try;f();catch ex;assert(strcmp(ex.identifier,'radia:cuda:ic'));return;end
error('radia:test:missingError','Invalid IC operation accepted');
end
