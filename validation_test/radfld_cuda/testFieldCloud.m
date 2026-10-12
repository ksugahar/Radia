function result=testFieldCloud()
% Lifecycle/error contracts, including stale tokens and an unloaded module.
x=[0 0 0;1 0 0];q=[1;-1];empty=zeros(0,1);
assert(radia_cuda_field('live_count')==0);
[handle,~]=radia_cuda_field('create_cloud',x,q,empty,empty,empty,1,0.05);
[a,~]=radia_cuda_field('field',handle,[2 1 0],'direct');
[b,~]=radia_cuda_field('field',handle,[2 1 0],'tree');
assert(norm(a-b)==0);
mustFail(@()checkedField(uint64(1),[2 1 0],'direct'));
mustFail(@()checkedField(double(handle),[2 1 0],'direct'));
mustFail(@()checkedField(handle,[NaN 1 0],'direct'));
mustFail(@()checkedField(handle,[2 1],'direct'));
mustFail(@()checkedField(handle,[2 1 0],'auto'));
radia_cuda_field('destroy',handle);
mustFail(@()radia_cuda_field('destroy',handle));
assert(radia_cuda_field('live_count')==0);
clear radia_cuda_field;
[newHandle,~]=radia_cuda_field('create_cloud',x,q,empty,empty,empty,1,0.05);
assert(newHandle~=handle);
mustFail(@()checkedField(handle,[2 1 0],'direct'));
[v,~]=radia_cuda_field('field',newHandle,zeros(0,3),'direct');assert(isequal(size(v),[0 3]));
radia_cuda_field('destroy',newHandle);
for i=1:5
    obj=radia.gpu.FieldCloud(x,q);obj.field([2 1 0],'direct');delete(obj);
end
assert(radia_cuda_field('live_count')==0);
result=struct('accepted',true,'live_handles',0);
end
function checkedField(handle,observations,algorithm)
[~,~]=radia_cuda_field('field',handle,observations,algorithm);
end
function mustFail(call)
failed=false;
try
    call();
catch exception
    failed=startsWith(exception.identifier,'radia:cuda:');
end
assert(failed,'Expected a checked native error');
end
