classdef FieldElements < handle
    %FIELDELEMENTS Persistent CUDA evaluator for retained analytic/curved RadFld elements.
    % Source is a version=1 struct with representation, volume and surface.
    % Coefficient records retain the CPU evaluator's conventions and ordering.
    % Curved sources also provide gauss_points and gauss_weights (column vectors).
    properties (SetAccess=private)
        SetupInfo
    end
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj = FieldElements(source,options)
            arguments
                source (1,1) struct
                options.LeafSize (1,1) double {mustBeInteger,mustBeInRange(options.LeafSize,1,65536)} = 32
                options.Theta (1,1) double {mustBePositive,mustBeLessThan(options.Theta,1)} = 0.05
                options.ImageMasks (:,1) double {mustBeInteger,mustBeInRange(options.ImageMasks,0,7)} = zeros(0,1)
                options.ImageSigns (:,1) double {mustBeReal,mustBeFinite} = zeros(0,1)
                options.ImageAngles (:,1) double {mustBeReal,mustBeFinite} = zeros(0,1)
            end
            timer=tic;
            if isempty(options.ImageAngles)
                options.ImageAngles=zeros(size(options.ImageMasks));
            end
            [obj.Token,obj.SetupInfo]=radia_cuda_field('create_elements',source, ...
                options.ImageMasks,options.ImageSigns,options.ImageAngles,options.LeafSize,options.Theta);
            obj.SetupInfo.total_seconds=toc(timer);
        end
        function [value,info] = field(obj,observations,algorithm)
            arguments
                obj (1,1) radia.gpu.FieldElements
                observations (:,3) double {mustBeReal,mustBeFinite}
                algorithm (1,:) char {mustBeMember(algorithm,{'direct','tree','tree_mixed32'})} = 'direct'
            end
            assert(obj.Token~=0,'radia:cuda:fieldDeleted','Field evaluator has been deleted');
            timer=tic;
            [value,info]=radia_cuda_field('field',obj.Token,observations,algorithm);
            info.total_seconds=toc(timer);
            info.algorithm=algorithm;
        end
        function delete(obj)
            if obj.Token~=0
                radia_cuda_field('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
end
