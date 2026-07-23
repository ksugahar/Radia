classdef GridFunction < handle
    %GRIDFUNCTION MATLAB owner for a native NGSolve GridFunction handle.
    properties (SetAccess=private)
        Handle (1,1) uint64 = uint64(0)
    end

    methods
        function obj = GridFunction(volPath, space, order, varargin)
            arguments
                volPath (1,1) string
                space (1,1) string
                order (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Repeating)
                varargin
            end
            obj.Handle = uint64(radia_mex('ngsolve.grid_function.create', ...
                char(volPath), char(space), order, varargin{:}));
        end

        function values = vector(obj, component)
            if nargin < 2, component = 1; end
            values = radia_mex('ngsolve.grid_function.vector', ...
                obj.Handle, component);
        end

        function coefficient = asCoefficient(obj)
            coefficient = radia.ngsolve.CoefficientFunction( ...
                radia_mex('ngsolve.grid_function.as_coefficient', obj.Handle), ...
                "native");
        end

        function delete(obj)
            if obj.Handle ~= 0
                radia_mex('ngsolve.grid_function.destroy', obj.Handle);
                obj.Handle = uint64(0);
            end
        end
    end

end
