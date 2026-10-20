function out=evaluateHDivBHChunks(chunks,law,m,wantTangent)
% Sum NGSolve-owned material chunks on GPU, gathering once at the solver boundary.
assert(iscell(chunks) && ~isempty(chunks));
validateattributes(wantTangent,{'logical'},{'scalar'});
nd=numel(m);assert(nd>0 && nd<=16384);
if ~isfield(law,'device')
    law.device=struct('knots',gpuArray(law.knots),'coefficients',gpuArray(law.coefficients),'integral_base',gpuArray(law.integral_base));
end
if isa(m,'gpuArray'),mDevice=m(:);else,mDevice=gpuArray(m(:));end
assert(gather(all(isfinite(mDevice))));
gradient=zeros(nd,1,'gpuArray');energy=gpuArray(0.);maximum=gpuArray(0.);
if wantTangent,tangent=gpuArray(sparse(nd,nd));end
for k=1:numel(chunks)
    q=chunks{k};assert(iscell(q.E) && numel(q.E)==3 && isa(q.weights,'gpuArray') && ~isempty(q.weights));
    for c=1:3,assert(isa(q.E{c},'gpuArray') && size(q.E{c},1)==numel(q.weights) && size(q.E{c},2)==nd);end
    partial=radia.gpu.evaluateHDivBHMaterial(q.E,q.weights,law,mDevice,wantTangent,true);
    gradient=gradient+partial.gradient;energy=energy+partial.energy;maximum=max(maximum,partial.maximum_magnetization);
    if wantTangent,tangent=tangent+partial.tangent;end
end
out=struct('gradient',gather(gradient),'energy',gather(energy),'maximum_magnetization',gather(maximum));
if wantTangent,out.tangent=gather(tangent);end
end
