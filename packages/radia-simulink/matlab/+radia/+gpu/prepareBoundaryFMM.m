function [operator,info] = prepareBoundaryFMM(volPath,boundaryLabels,options)
%PREPAREBOUNDARYFMM Native flat-P1 Galerkin maps with exact near corrections.
% The caller supplies a check-vol-approved exterior .vol. NGSolve owns P1
% basis, normals and quadrature. Shared Sauter-Schwab integration corrects
% all near pairs before the FMM approximation is applied to the remaining
% point-rule interactions. Unsupported geometry or interfaces raise.
arguments
    volPath (1,1) string
    boundaryLabels (1,:) string
    options.Kind (1,1) string {mustBeMember(options.Kind,["single-layer","double-layer"])} = "single-layer"
    options.Order (1,1) double {mustBeInteger,mustBeInRange(options.Order,3,8)} = 8
    options.Separation (1,1) double {mustBeInRange(options.Separation,1.1,8)} = 2
end
timer=tic;maps=radia_cuda_fmm('prepare_boundary',char(volPath),cellstr(boundaryLabels));
n=maps.dofs;ns=size(maps.points,1);
source=gpuArray(sparse(maps.map_rows+1,maps.map_columns+1,maps.map_values,ns,n));
if options.Kind=="single-layer";values=maps.correction_sl;else;values=maps.correction_dl;end
correction=gpuArray(sparse(maps.correction_rows+1,maps.correction_columns+1,values,n,n));
operator=radia.gpu.LaplaceBoundaryFMM(maps.points,maps.points,source,source',correction, ...
    Kind=options.Kind,Normals=maps.normals,Order=options.Order,Separation=options.Separation,SelfInteraction="omit");
info=struct('schema',maps.schema,'global_dofs',maps.global_dofs,'boundary_dofs',n,'quadrature_points',ns, ...
    'exact_near_pairs',maps.near_pairs,'setup_seconds',toc(timer),'near_quadrature','Shared P1 Sauter-Schwab, singular order 10; regular degree 17', ...
    'far_quadrature','NGSolve triangle degree 16','geometry','Flat exterior scalar P1 only');
end
