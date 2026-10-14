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
% Exercise element handle creation/destruction as well as malformed schemas.
% Independent numerical parity is checked separately by validate.py.
vertices=[0 0 0;1 0 0;0 1 0;0 0 1];
kinds={'tet','polynomial_tet','cubic_tet','sextic_quartic_tet','curved_tet'};
widths=[16 22;32 22;32 29;96 44;34 24];
for k=1:numel(kinds)
    source=struct('version',1,'representation',kinds{k}, ...
        'volume',zeros(1,widths(k,1)),'surface',zeros(0,widths(k,2)));
    if strcmp(kinds{k},'curved_tet')
        edges=[1 2;2 3;3 1;1 4;2 4;3 4];
        nodes=[vertices;(vertices(edges(:,1),:)+vertices(edges(:,2),:))/2];
        source.volume(1:30)=reshape(nodes.',1,30);
        source.volume(31)=1;
        source.gauss_points=[0.211324865405187;0.788675134594813];
        source.gauss_weights=[0.5;0.5];
    else
        source.volume(1:12)=reshape(vertices.',1,12);source.volume(13)=1;
    end
    object=radia.gpu.FieldElements(source,LeafSize=1);
    [value,~]=object.field([2 3 4],'direct');assert(all(isfinite(value),'all'));
    delete(object);
    malformed=source;malformed.version=2;
    mustFail(@()checkedElements(malformed));
    malformed=source;malformed.volume(1)=NaN;
    mustFail(@()checkedElements(malformed));
    malformed=source;malformed.surface=zeros(0,1);
    mustFail(@()checkedElements(malformed));
    if strcmp(kinds{k},'curved_tet')
        malformed=source;malformed.gauss_weights=[-1;2];
        mustFail(@()checkedElements(malformed));
    end
    assert(radia_cuda_field('live_count')==0);
end
result=struct('accepted',true,'live_handles',0,'element_families',{kinds});
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

function checkedElements(source)
empty=zeros(0,1);
[~,~]=radia_cuda_field('create_elements',source,empty,empty,empty,1,0.05);
end
