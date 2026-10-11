function [solution,info]=amsPCG(snapshot,rhs,tolerance,maxIterations,applyOnly)
%AMSPCG Experimental FP64 assembled AMS on caller-owned GPU arrays.
% CPU-owned NGSolve setup is exchanged through readAMSExchange. Upload the
% numeric hierarchy explicitly; this entry point never substitutes a CPU solve.
if nargin<3;tolerance=1e-9;end
if nargin<4;maxIterations=2000;end
if nargin<5;applyOnly=false;end
assert(exist('radia_cuda_ams','file')==3,'radia:gpu:missingAMS', ...
    'AMS CUDA MEX is missing; build/install the matching radia-gpu artifact');
[solution,info]=radia_cuda_ams(snapshot,rhs,tolerance,maxIterations,applyOnly);
end
