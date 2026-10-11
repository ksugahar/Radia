function result = check_gateway(binaryDir)
% Native boundary checks independent of finite-element Poisson validation.
addpath(binaryDir);
rp=gpuArray(int32([0;2;5;7])); ci=gpuArray(int32([0;1;0;1;2;1;2]));
av=gpuArray([2;-1;-1;2;-1;-1;2]); b=gpuArray([1;0;1]); di=gpuArray([.5;.5;.5]);
[x,info]=radia_cuda_csr_pcg(rp,ci,av,b,di,1e-12,100);
assert(norm(gather(x)-ones(3,1))<1e-12);
assert(info.true_relative_residual<1e-12);
assert(rejected(@() twoOutputs(gather(rp),ci,av,b,di,1e-12,100)));
assert(rejected(@() twoOutputs(rp,gpuArray(int32([0;1;0;3;2;1;2])),av,b,di,1e-12,100)));
assert(rejected(@() twoOutputs(rp,ci,av,b,gpuArray([.5;0;.5]),1e-12,100)));
assert(rejected(@() twoOutputs(rp,ci,av,b,di,NaN,100)));
assert(rejected(@() twoOutputs(rp,ci,av,b,di,1e-12,1)));
assert(rejected(@() twoOutputs(rp,ci,av,gpuArray([1;NaN;1]),di,1e-12,100)));
[zero,info]=radia_cuda_csr_pcg(rp,ci,av,gpuArray(zeros(3,1)),di,1e-12,100);
assert(norm(gather(zero))==0 && info.iterations==0);
for k=1:10
    [x,~]=radia_cuda_csr_pcg(rp,ci,av,b,di,1e-12,100);
    assert(norm(gather(x)-ones(3,1))<1e-12);
end
clear radia_cuda_csr_pcg;
[x,~]=radia_cuda_csr_pcg(rp,ci,av,b,di,1e-12,100);
assert(norm(gather(x)-ones(3,1))<1e-12);
result=struct('passed',true,'cases',10,'repeat_solves',10,'clear_and_reload',true);
end
function twoOutputs(varargin)
[~,~]=radia_cuda_csr_pcg(varargin{:});
end
function ok=rejected(action)
try
    action(); ok=false;
catch exception
    ok=startsWith(exception.identifier,'radia:cuda:');
end
end
