function [model,fem,info] = assembleHCurl(volPath,boundaryLabels,options)
%ASSEMBLEHCURL Build affine tet6 magnetic operators within MATLAB/mexcuda.
% NGSolve owns orientation, geometry, quadrature, shapes and gradient.
% CUDA assembles nu*curl(u)*curl(v)+MassCoefficient*u*v and J*v in SI.
% Positive mass fixes the gauge. This is a bounded formulation, not a general
% NGSolve symbolic-form adapter. Curved/high-order/non-tet elements raise.
% model is the native prepareAMS input; fem is GPU-resident nonlinear input.
arguments
    volPath
    boundaryLabels (1,:) string
    options.Dirichlet (1,1) string = ".*"
    options.Reluctivity (1,1) double {mustBeFinite,mustBePositive} = 1
    options.MassCoefficient (1,1) double {mustBeFinite,mustBePositive} = 1
    options.CurrentDensity (1,3) double {mustBeFinite} = [1 0 0]
    options.Gamma (1,1) double {mustBeFinite,mustBeNonnegative} = 0
end
if isa(volPath,'radia.gpu.HCurlGeometry')
    assert(isequal(sort(boundaryLabels),sort(volPath.BoundaryLabels)) && options.Dirichlet==volPath.Dirichlet, ...
        'radia:cuda:geometryContract','Geometry snapshot labels or constraints differ');
    [model,fem,info]=volPath.assemble(Reluctivity=options.Reluctivity, ...
        MassCoefficient=options.MassCoefficient,CurrentDensity=options.CurrentDensity,Gamma=options.Gamma);
    return;
end
assert((isstring(volPath)&&isscalar(volPath)) || (ischar(volPath)&&isrow(volPath)), ...
    'radia:cuda:geometryType','Expected a .vol path or HCurlGeometry object');
assert(isfile(volPath),'radia:cuda:vol','Checked .vol path does not exist');
[model,fem,info]=radia_cuda_ams('assemble_hcurl',char(volPath),cellstr(boundaryLabels), ...
    char(options.Dirichlet),options.Reluctivity,options.MassCoefficient,options.CurrentDensity,options.Gamma);
end
