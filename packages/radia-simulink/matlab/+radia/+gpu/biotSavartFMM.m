function [B,A,info] = biotSavartFMM(context,weightedCurrent,options)
%BIOTSAVARTFMM Free-space MQS field from explicitly weighted current samples.
% weightedCurrent is J*dV or I*dl in SI, ns x 3 gpuArray double.
% Quadrature/filament near corrections remain caller-owned. No wire-radius,
% sampling rule or coincident self field is inferred by this function.
arguments
    context (1,1) radia.gpu.LaplaceFMM
    weightedCurrent
    options.Permeability (1,1) double {mustBeFinite,mustBePositive} = 1.25663706212e-6
end
assert(isa(weightedCurrent,'gpuArray')&&size(weightedCurrent,2)==3, ...
    'radia:fmm:current','Weighted current must be an ns x 3 gpuArray');
[potential,field,info]=context.apply(weightedCurrent);
scale=options.Permeability/(4*pi);
A=scale*potential;
B=scale*[field(:,6)-field(:,8),field(:,7)-field(:,3),field(:,2)-field(:,4)];
info.physical_scaling='SI: permeability/(4*pi); A in T*m, B in T';
end
