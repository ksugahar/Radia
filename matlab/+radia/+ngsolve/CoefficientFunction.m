classdef CoefficientFunction < handle
    %COEFFICIENTFUNCTION MATLAB owner for a native NGSolve coefficient handle.
    properties (SetAccess=private)
        Handle (1,1) uint64 = uint64(0)
    end

    methods
        function obj = CoefficientFunction(values, mode)
            if nargin == 2 && mode == "native"
                obj.Handle = uint64(values);
            else
                obj.Handle = uint64(radia_mex( ...
                    'ngsolve.coefficient_function.constant_create', values));
            end
        end

        function info = info(obj)
            info = radia_mex('ngsolve.coefficient_function.info', obj.Handle);
        end

        function delete(obj)
            if obj.Handle ~= 0
                radia_mex('ngsolve.coefficient_function.destroy', obj.Handle);
                obj.Handle = uint64(0);
            end
        end
    end

end
