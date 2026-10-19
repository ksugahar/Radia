classdef LaplaceFMM < handle
    %LAPLACEFMM Numeric free-space 1/r scalar/dipole CUDA FMM.
    % All geometry/tree setup and execution complete within MATLAB/mexcuda.
    % Potential has no 1/(4*pi) factor; Field is -grad(Potential).
    % charge is ns x channels, dipole is ns x (3*channels), channels=1..3.
    % Coincident points raise unless SelfInteraction="omit" is explicit.
    % Surface quadrature, exact near corrections and physical scaling belong
    % to the calling operator. This class does not infer a BEM discretization.
    properties (SetAccess=private,Hidden)
        Token uint64 = uint64(0)
        SourceCount double
    end
    methods
        function obj=LaplaceFMM(sources,targets,options)
            arguments
                sources (:,3) double {mustBeFinite}
                targets (:,3) double {mustBeFinite}
                options.Order (1,1) double {mustBeInteger,mustBeInRange(options.Order,3,8)} = 6
                options.LeafCapacity (1,1) double {mustBeInteger,mustBeInRange(options.LeafCapacity,8,512)} = 64
                options.Separation (1,1) double {mustBeInRange(options.Separation,1.1,8)} = 1.5
                options.SelfInteraction (1,1) string {mustBeMember(options.SelfInteraction,["error","omit"])} = "error"
            end
            obj.SourceCount=size(sources,1);
            obj.Token=radia_cuda_fmm('create',sources,targets,options.Order,options.LeafCapacity, ...
                options.Separation,options.SelfInteraction=="omit");
        end
        function delete(obj)
            if obj.Token~=0
                token=obj.Token;obj.Token=uint64(0);radia_cuda_fmm('destroy',token);
            end
        end
        function [potential,field,info]=apply(obj,charge,dipole)
            if nargin<3;dipole=zeros(obj.SourceCount,3*size(charge,2),'like',charge);end
            assert(isa(charge,'gpuArray')&&isa(dipole,'gpuArray'), ...
                'radia:fmm:backend','Densities must be explicit gpuArray values');
            assert(size(charge,1)==obj.SourceCount&&size(dipole,1)==obj.SourceCount, ...
                'radia:fmm:shape','Density rows must match source count');
            [potential,field,info]=radia_cuda_fmm('apply',obj.Token,charge,dipole);
        end
    end
end
