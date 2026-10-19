function [field,info] = fieldCloudFMM(context,strength)
%FIELDCLOUDFMM Explicit FMM route for RadFld retained numeric cloud sources.
% sum q*(target-source)/r^3, without 1/(4*pi). The context owns fixed targets.
% Analytic/curved element sources and image transforms are not sampled here.
arguments
    context (1,1) radia.gpu.LaplaceFMM
    strength
end
assert(isa(strength,'gpuArray')&&size(strength,2)==1, ...
    'radia:fmm:strength','Cloud strength must be an ns x 1 gpuArray');
[~,field,info]=context.apply(strength);
info.convention='Existing RadFld retained-cloud convention, without 1/(4*pi)';
end
