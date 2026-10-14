function fmm_boundary(volPath,binaryPath,matlabPath,outputPath)
addpath(binaryPath);addpath(matlabPath);initial=radia_cuda_fmm('live_count');
labels=["back","front","left","right","top","bottom"];
[sl,info]=radia.gpu.prepareBoundaryFMM(volPath,labels);slCleanup=onCleanup(@()delete(sl));
[dl,dlInfo]=radia.gpu.prepareBoundaryFMM(volPath,labels,Kind="double-layer");dlCleanup=onCleanup(@()delete(dl));
n=info.boundary_dofs;x=sin((1:n)'*.37)+.1*cos((1:n)'*.13);density=gpuArray(x);
[v,slApply]=sl.apply(density);slResult=gather(v);[v,dlApply]=dl.apply(density);dlResult=gather(v);
global_dofs=info.global_dofs;report=struct('single_layer',info,'double_layer',dlInfo,'sl_apply',slApply,'dl_apply',dlApply);
clear slCleanup dlCleanup sl dl;assert(radia_cuda_fmm('live_count')==initial);
save(outputPath,'x','slResult','dlResult','global_dofs','report');
end
