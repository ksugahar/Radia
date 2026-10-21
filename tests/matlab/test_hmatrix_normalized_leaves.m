function [result,fixture]=test_hmatrix_normalized_leaves()
% Explicit normalized ABI: flat factors, one-based offsets/permutation.
% Dense values are row-major; low-rank factors are column-major [size,rank].
gpuDevice();
before=radia_cuda_hmatrix('live_count');
leaf=struct('row',1,'column',1,'rows',2,'columns',2,'rank',0,'a1',[3;1;1;4],'a2',[]);
leaves=repmat(leaf,1,3);
leaves(2)=struct('row',1,'column',3,'rows',2,'columns',1,'rank',1,'a1',2,'a2',[.5;-.25]);
leaves(3)=struct('row',3,'column',3,'rows',1,'columns',1,'rank',0,'a1',5,'a2',[]);
cfg=struct('schema_version',2,'assembly','normalized_leaves','dimension',3, ...
    'permutation',[3;1;2],'scale',[2;3;4],'leaves',leaves,'operator_kind','gram');
% Independent explicit normalized matrix, no normalization or dense reconstruction of factors.
G=[3,1,1;1,4,-.5;1,-.5,5]; p=cfg.permutation; scale=cfg.scale;
x=[.25;-2;3]; gathered=x(p).*scale(p); expected=zeros(3,1); expected(p)=(G*gathered).*scale(p);
op=radia.gpu.HMatrixOperator(cfg); op.upload(x); op.apply(); actual=op.download(); delete(op);
err=norm(actual-expected)/norm(expected); assert(err<1e-14);
% Unknown schema/options, duplicate ordering, lower leaves and malformed factors fail without handles.
bad=cfg; bad.schema_version=1; reject(bad);
bad=cfg; bad.extra=true; reject(bad);
bad=cfg; bad.permutation=[1;1;3]; reject(bad);
bad=cfg; bad.leaves(2).row=3; bad.leaves(2).column=1; reject(bad);
bad=cfg; bad.leaves(2).a2=1; reject(bad);
bad=cfg; bad.leaves=[leaves,leaves(2)]; reject(bad); % exact duplicate
bad=cfg; bad.leaves=leaves([1,3]); reject(bad); % missing off-diagonal block
bad=cfg; bad.leaves=[leaves,struct('row',1,'column',2,'rows',1,'columns',1,'rank',0,'a1',1,'a2',[])]; reject(bad); % partial overlap
op=radia.gpu.HMatrixOperator(cfg); op.upload(x); op.apply(); assert(norm(op.download()-expected)<1e-12); delete(op);
assert(radia_cuda_hmatrix('live_count')==before);
result=struct('accepted',true,'relative_error',err,'live_count',before);
fixture=struct('config',cfg,'x',x,'expected',expected);
end
function reject(cfg)
failed=false; try, op=radia.gpu.HMatrixOperator(cfg); delete(op); catch, failed=true; end
assert(failed,'Expected invalid normalized payload rejection');
end
