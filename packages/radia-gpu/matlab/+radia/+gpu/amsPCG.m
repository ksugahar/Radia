function [solution,info]=amsPCG(snapshot,rhs,tolerance,maxIterations,applyOnly,useGraph,currentOperator)
%AMSPCG Experimental FP64 assembled AMS on caller-owned GPU arrays.
% Optional currentOperator is a same-size SPD int32 CSR GPU tangent with
% constrained rows/columns removed. snapshot remains the frozen preconditioner.
% PCG and its true residual use currentOperator; no implicit rebuild/fallback.
% MATLAB/mexcuda setup is available through prepareAMS. Upload the
% numeric hierarchy explicitly; this entry point never substitutes a CPU solve.
% applyOnly checks finite output but does not solve a linear system:
% info.residual_checked is false and info.true_relative_residual is NaN.
if nargin<3;tolerance=1e-9;end
if nargin<4;maxIterations=2000;end
if nargin<5;applyOnly=false;end
% Experimental Graph replay is opt-in until end-to-end benefit is established.
if nargin<6;useGraph=false;end
assert(exist('radia_cuda_ams','file')==3,'radia:gpu:missingAMS', ...
    'AMS CUDA MEX is missing; build/install the matching radia-gpu artifact');
if nargin>=7
    [solution,info]=radia_cuda_ams(snapshot,rhs,tolerance,maxIterations,applyOnly,useGraph,currentOperator);
else
    [solution,info]=radia_cuda_ams(snapshot,rhs,tolerance,maxIterations,applyOnly,useGraph);
end
end
