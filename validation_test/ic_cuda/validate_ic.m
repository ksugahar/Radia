function validate_ic(inputFile,binaryPath,matlabPath,outputFile)
% Independent Python numerical fixture, native rejection and ownership gates.
addpath(matlabPath);addpath(binaryPath,'-begin');data=load(inputFile);
before=radia_cuda_ic('live_count');op=toCSR(data.A);rhs=gpuArray(data.b);
[h,setup]=radia_cuda_ic('create',op,1.05);cleanup=onCleanup(@()radia_cuda_ic('release',h));
applied=gather(radia_cuda_ic('apply',h,rhs));applyError=norm(applied-data.expected_pre)/norm(data.expected_pre);
assert(isfinite(applyError)&&applyError<1e-12);
scales=[1 2 .5 0 1];errors=zeros(5,2);reports=cell(5,1);
for k=1:5
    [x,reports{k}]=radia_cuda_ic('solve',h,scales(k)*rhs,1e-12,100);
    x=gather(x);errors(k,:)=[norm(x-scales(k)*data.expected)/max(norm(scales(k)*data.expected),1), ...
        norm(data.A*x-scales(k)*data.b)/max(norm(scales(k)*data.b),1)];
end
assert(all(isfinite(errors),'all')&&max(errors,[],'all')<1e-11);
sparseApplyError=[];
if isfield(data,'sparse_A')
    sparseOp=toSparseCSR(data.sparse_A);[sh,~]=radia_cuda_ic('create',sparseOp,1.05);
    sparseCleanup=onCleanup(@()radia_cuda_ic('release',sh));
    pre=gather(radia_cuda_ic('apply',sh,gpuArray(data.sparse_b)));
    sparseApplyError=norm(pre-data.sparse_expected_pre)/norm(data.sparse_expected_pre);
    [sx,~]=radia_cuda_ic('solve',sh,gpuArray(data.sparse_b),1e-12,1000);sx=gather(sx);
    assert(isfinite(sparseApplyError)&&sparseApplyError<1e-11);
    assert(norm(sx-data.sparse_expected)/norm(data.sparse_expected)<1e-10);
    assert(norm(data.sparse_A*sx-data.sparse_b)/norm(data.sparse_b)<1e-10);clear sparseCleanup;
end
object=radia.gpu.ICContext(op,Shift=1.05);objectCleanup=onCleanup(@()delete(object));
[objectSolution,objectInfo]=object.solve(rhs,Tolerance=1e-12,ReturnOnGPU=false);
assert(norm(objectSolution-data.expected)<1e-11&&objectInfo.result_download_bytes==24);
assert(norm(gather(object.apply(rhs))-data.expected_pre)<1e-11);clear objectCleanup;
rejections=0;
masked=toCSR(blkdiag(data.A,0));masked.free=gpuArray(int32([1;1;1;0]));
[mh,maskInfo]=radia_cuda_ic('create',masked,1.05);maskCleanup=onCleanup(@()radia_cuda_ic('release',mh));
[mx,~]=radia_cuda_ic('solve',mh,gpuArray([data.b;0]),1e-12,100);
assert(norm(gather(mx)-[data.expected;0])<1e-11&&maskInfo.explicit_free_mask);
assert(norm(gather(masked.values)-reshape(blkdiag(data.A,0).',[],1))==0);
mustFail(@()solveLimited(mh,gpuArray([data.b;1])));rejections=rejections+1;clear maskCleanup;
bad=op;bad.free=gpuArray(int32([1;2;1]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
bad=op;bad.free=gpuArray(int32([1;1;0]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
mustFail(@()solveLimited(h,rhs));rejections=rejections+1;
mustFail(@()solveLimited(uint64(1),rhs));rejections=rejections+1;
mustFail(@()solveLimited(double(h),rhs));rejections=rejections+1;
mustFail(@()solveLimited(h,gpuArray([NaN;1;1])));rejections=rejections+1;
mustFail(@()solveLimited(h,gpuArray([Inf;1;1])));rejections=rejections+1;
mustFail(@()solveExtreme(h,gpuArray([1e-300;0;0])));rejections=rejections+1;
tiny=realmin('double')*eps('double');assert(tiny~=0);
mustFail(@()solveExtreme(h,gpuArray([tiny;0;0])));rejections=rejections+1;
% A rejected solve does not invalidate a live context.
[recovered,~]=radia_cuda_ic('solve',h,rhs,1e-12,100);assert(norm(gather(recovered)-data.expected)<1e-11);
for A={ [1 2;2 1], [1 1;2 3], [NaN 0;0 1], [Inf 0;0 1], [-1 0;0 1] }
    mustFail(@()badCreate(toCSR(A{1}),1));rejections=rejections+1;
end
bad=op;bad.columns=gpuArray(int32([0;0;2;0;1;2;0;1;2]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
bad=op;bad.columns=gpuArray(int32([0;1;3;0;1;2;0;1;2]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
bad=op;bad.rowptr=gpuArray(int32([1;3;6;9]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
bad=op;bad.values=gather(bad.values);mustFail(@()badCreate(bad,1));rejections=rejections+1;
bad=struct('rows',2.,'cols',2.,'rowptr',gpuArray(int32([0;1;2])), ...
    'columns',gpuArray(int32([1;0])),'values',gpuArray([1.;1.]));mustFail(@()badCreate(bad,1));rejections=rejections+1;
mustFail(@()badCreate(op,NaN));rejections=rejections+1;
mustFail(@()badCreate(op,.5));rejections=rejections+1;
clear cleanup;
mustFail(@()solveLimited(h,rhs));rejections=rejections+1;
mustFail(@()radia_cuda_ic('release',h));rejections=rejections+1;
% Generations distinguish a released object from a fresh one.
[fresh,~]=radia_cuda_ic('create',op,1.05);assert(fresh~=h);radia_cuda_ic('release',fresh);
assert(radia_cuda_ic('live_count')==before);
result=struct('passed',true,'apply_relative',applyError,'sparse_apply_relative',sparseApplyError,'solve_errors',errors, ...
    'rejection_count',rejections,'setup',setup,'solve',{reports},'live_count_restored',true, ...
    'loaded_binary',which('radia_cuda_ic'),'matlab',version,'device',gpuDevice().Name);
save(outputFile,'result','-v7');
end
function op=toCSR(A)
% Deliberately retain the full dense pattern, including structural zeros.
n=size(A,1);op=struct('rows',double(n),'cols',double(n), ...
    'rowptr',gpuArray(int32((0:n:n*n)')),'columns',gpuArray(repmat(int32((0:n-1)'),n,1)), ...
    'values',gpuArray(reshape(A.',[],1)));
end
function solveLimited(h,b)
[~,~]=radia_cuda_ic('solve',h,b,1e-12,1);
end
function op=toSparseCSR(A)
n=size(A,1);[columns,rows,values]=find(sparse(A.'));counts=accumarray(rows,1,[n 1]);
op=struct('rows',double(n),'cols',double(n),'rowptr',gpuArray(int32([0;cumsum(counts)])), ...
    'columns',gpuArray(int32(columns-1)),'values',gpuArray(values));
end
function badCreate(op,shift)
[h,~]=radia_cuda_ic('create',op,shift);radia_cuda_ic('release',h);
end
function solveExtreme(h,b)
[~,~]=radia_cuda_ic('solve',h,b,1e-10,100);
end
function mustFail(callback)
failed=false;try;callback();catch ex;assert(startsWith(ex.identifier,'radia:cuda:'),ex.identifier);failed=true;end
assert(failed,'radia:test:missingError','Invalid IC operation was accepted');
end
