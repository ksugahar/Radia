function fmm_physics(inputPath,binaryPath,matlabPath,outputPath)
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);initial=radia_cuda_fmm('live_count');
map=gpuArray(spdiags(s.weights,0,numel(s.weights),numel(s.weights)));
test=gpuArray(speye(size(s.targets,1)));correction=gpuArray(sparse(size(s.targets,1),numel(s.weights)));
sl=radia.gpu.LaplaceBoundaryFMM(s.sphere,s.targets,map,test,correction,Order=8,LeafCapacity=32,Separation=2);
slCleanup=onCleanup(@()delete(sl));
[v,slInfo]=sl.apply(gpuArray(ones(numel(s.weights),1)));slPotential=gather(v);clear slCleanup sl;
dl=radia.gpu.LaplaceBoundaryFMM(s.sphere,s.targets,map,test,correction,Kind="double-layer", ...
    Normals=s.sphere,Order=8,LeafCapacity=32,Separation=2);
dlCleanup=onCleanup(@()delete(dl));
[v,dlInfo]=dl.apply(gpuArray(ones(numel(s.weights),1)));dlPotential=gather(v);clear dlCleanup dl;
context=radia.gpu.LaplaceFMM(s.coil,s.coilTargets,Order=8,LeafCapacity=32,Separation=2);
contextCleanup=onCleanup(@()delete(context));
[b,a,bsInfo]=radia.gpu.biotSavartFMM(context,gpuArray(s.current));B=gather(b);A=gather(a);
clear contextCleanup context;
cloud=radia.gpu.LaplaceFMM(s.sphere,s.targets,Order=8,LeafCapacity=32,Separation=2);
cloudCleanup=onCleanup(@()delete(cloud));
[v,cloudInfo]=radia.gpu.fieldCloudFMM(cloud,gpuArray(s.weights));cloudField=gather(v);
% Failed density validation must leave the context reusable.
bad=gpuArray(s.weights);bad(1)=NaN;caught=false;
try;[~,~,~]=cloud.apply(bad);catch error;caught=startsWith(error.identifier,'radia:cuda:fmm');end
assert(caught,'radia:fmm:nan','Nonfinite density was not rejected');
[v,~]=radia.gpu.fieldCloudFMM(cloud,gpuArray(s.weights));assert(isequal(cloudField,gather(v)));
token=cloud.Token;clear cloudCleanup cloud;caught=false;
try;radia_cuda_fmm('destroy',token);catch error;caught=startsWith(error.identifier,'radia:cuda:fmm');end
assert(caught&&radia_cuda_fmm('live_count')==initial,'radia:fmm:lifecycle','Stale handle or context leak');
% Coincident point policy is explicit; failed apply remains destructible.
self=radia.gpu.LaplaceFMM(s.sphere,s.sphere(1:2,:),Order=4);selfCleanup=onCleanup(@()delete(self));caught=false;
try;[~,~,~]=self.apply(gpuArray(s.weights));catch error;caught=contains(error.message,'Coincident');end
assert(caught,'radia:fmm:self','Undeclared self interaction was not rejected');clear selfCleanup self;
assert(radia_cuda_fmm('live_count')==initial);
% This validation owns its Engine and must not unload caller-owned contexts.
assert(initial==0,'radia:fmm:isolated','Reload validation requires an isolated registry');
old=radia_cuda_fmm('create',s.sphere,s.targets,8,32,2,false);
radia_cuda_fmm('destroy',old);clear radia_cuda_fmm;
fresh=radia_cuda_fmm('create',s.sphere,s.targets,8,32,2,false);
freshCleanup=onCleanup(@()radia_cuda_fmm('destroy',fresh));
assert(old~=fresh,'radia:fmm:generation','Reload recycled an old token');
caught=false;
try;[~,~,~]=radia_cuda_fmm('apply',old,gpuArray(s.weights),gpuArray(zeros(numel(s.weights),3)));
catch error;caught=startsWith(error.identifier,'radia:cuda:fmm');end
assert(caught,'radia:fmm:reload','Old apply handle survived reload');
caught=false;
try;radia_cuda_fmm('destroy',old);catch error;caught=startsWith(error.identifier,'radia:cuda:fmm');end
assert(caught,'radia:fmm:reload','Old destroy handle survived reload');
[v,e,~]=radia_cuda_fmm('apply',fresh,gpuArray(s.weights),gpuArray(zeros(numel(s.weights),3)));
assert(all(isfinite(gather(v)),'all')&&all(isfinite(gather(e)),'all'));
clear freshCleanup;
assert(radia_cuda_fmm('live_count')==0);
report=struct('sl',slInfo,'dl',dlInfo,'biot_savart',bsInfo,'cloud',cloudInfo, ...
    'nan_rejected',true,'stale_handle_rejected',true,'reload_stale_handle_rejected',true, ...
    'self_rejected',true,'live_count_restored',true);
save(outputPath,'slPotential','dlPotential','B','A','cloudField','report');
end
