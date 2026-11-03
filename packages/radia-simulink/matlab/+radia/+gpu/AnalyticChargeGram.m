classdef AnalyticChargeGram < handle
    %ANALYTICCHARGEGRAM Device-filled Gram of constant-density tetrahedra and triangles.
    % cells and faces are flat row-major coordinate columns: 12 per linear
    % tetrahedron, 9 per linear triangle, at most 128 charges. The block is
    % filled on the device from this geometry with the all-analytic near route;
    % only vectors cross. No images, curved or high-order geometry, charge map
    % or material solve. It is a bounded block, not an HDiv-MMM solver.
    properties (Access=private)
        Token = uint64(0)
    end
    methods
        function obj=AnalyticChargeGram(cells,faces)
            arguments
                cells (:,1) double {mustBeReal,mustBeFinite}
                faces (:,1) double {mustBeReal,mustBeFinite}
            end
            obj.Token=radia_hdiv_device('create_analytic_gram',cells,faces);
        end
        function y=apply(obj,x)
            arguments
                obj (1,1) radia.gpu.AnalyticChargeGram
                x (:,1) double {mustBeReal,mustBeFinite}
            end
            y=radia_hdiv_device('apply',obj.Token,x);
        end
        function info=assembly(obj)
            info=radia_hdiv_device('assembly',obj.Token);
        end
        function info=stats(obj)
            info=radia_hdiv_device('stats',obj.Token);
        end
        function delete(obj)
            if obj.Token~=0
                radia_hdiv_device('destroy',obj.Token);
                obj.Token=uint64(0);
            end
        end
    end
end
