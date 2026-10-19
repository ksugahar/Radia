classdef LaplaceBoundaryFMM < handle
    %LAPLACEBOUNDARYFMM Explicit Galerkin maps around numeric Laplace FMM.
    % y=TargetMap*potential(SourceMap*x)/(4*pi)+NearCorrection*x.
    % SourceMap includes quadrature weights and trial shape evaluation;
    % TargetMap includes test shapes and target quadrature weights.
    % NearCorrection replaces the point-rule near/self blocks with checked
    % analytic/singular quadrature blocks. It is required, never inferred.
    % These maps may compose the boundary block of a FEM-BEM system.
    properties (Access=private)
        Context
        SourceMap
        TargetMap
        NearCorrection
        Normals
        Kind
    end
    methods
        function obj=LaplaceBoundaryFMM(sources,targets,sourceMap,targetMap,nearCorrection,options)
            arguments
                sources (:,3) double {mustBeFinite}
                targets (:,3) double {mustBeFinite}
                sourceMap
                targetMap
                nearCorrection
                options.Kind (1,1) string {mustBeMember(options.Kind,["single-layer","double-layer"])} = "single-layer"
                options.Normals (:,3) double {mustBeFinite} = zeros(0,3)
                options.Order (1,1) double {mustBeInteger,mustBeInRange(options.Order,3,8)} = 6
                options.LeafCapacity (1,1) double {mustBeInteger,mustBeInRange(options.LeafCapacity,8,512)} = 64
                options.Separation (1,1) double {mustBeInRange(options.Separation,1.1,8)} = 1.5
                options.SelfInteraction (1,1) string {mustBeMember(options.SelfInteraction,["error","omit"])} = "error"
            end
            assert(isa(sourceMap,'gpuArray')&&isa(targetMap,'gpuArray')&&isa(nearCorrection,'gpuArray'), ...
                'radia:fmm:maps','Explicit GPU source, target and near-correction maps are required');
            assert(size(sourceMap,1)==size(sources,1)&&size(targetMap,2)==size(targets,1)&& ...
                isequal(size(nearCorrection),[size(targetMap,1),size(sourceMap,2)]), ...
                'radia:fmm:maps','Incompatible boundary map dimensions');
            if options.Kind=="double-layer"
                assert(size(options.Normals,1)==size(sources,1)&&all(abs(vecnorm(options.Normals,2,2)-1)<1e-12), ...
                    'radia:fmm:normals','Double layer requires explicit oriented unit source normals');
            end
            obj.SourceMap=sourceMap;obj.TargetMap=targetMap;obj.NearCorrection=nearCorrection;
            obj.Normals=gpuArray(options.Normals);obj.Kind=options.Kind;
            obj.Context=radia.gpu.LaplaceFMM(sources,targets,Order=options.Order, ...
                LeafCapacity=options.LeafCapacity,Separation=options.Separation,SelfInteraction=options.SelfInteraction);
        end
        function delete(obj)
            if ~isempty(obj.Context);delete(obj.Context);obj.Context=[];end
        end
        function [value,info]=apply(obj,x)
            assert(isa(x,'gpuArray')&&size(x,1)==size(obj.SourceMap,2)&&size(x,2)==1, ...
                'radia:fmm:density','Boundary density must be a checked GPU coefficient vector');
            density=obj.SourceMap*x;
            if obj.Kind=="single-layer"
                [potential,~,info]=obj.Context.apply(density);
            else
                dipole=density.*obj.Normals;
                [potential,~,info]=obj.Context.apply(zeros(size(density),'like',density),dipole);
            end
            value=(obj.TargetMap*potential)/(4*pi)+obj.NearCorrection*x;
            assert(gather(all(isfinite(value))),'radia:fmm:finite','Nonfinite corrected boundary result');
        end
    end
end
