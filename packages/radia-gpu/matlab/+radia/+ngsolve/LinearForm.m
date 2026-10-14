classdef LinearForm < handle
    %LINEARFORM Persistent native NGSolve LinearForm and RHS vector.
    %   The MEX contract provides constant volume sources and real
    %   CoefficientFunction volume sources for H1, HCurl, and HDiv spaces.

    properties (SetAccess=private)
        Space string = ""
        Source string = ""
        Label string = ""
        Size double = 0
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = LinearForm(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Space = string(info.space);
            obj.Source = string(info.source);
            obj.Label = string(info.label);
            obj.Size = info.size;
        end
    end

    methods (Static)
        function obj = create(space, source, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                source (1,1) string = "constant"
                options.Value (1,1) double = 1.0
                options.Label (1,1) string = "radia_matlab_linear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.linear_form.create', space.nativeHandle(), ...
                char(source), options.Value, char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.linear_form.info', nativeHandle);
            obj = radia.ngsolve.LinearForm(nativeHandle, info);
        end

        function obj = createFromCoefficient(space, coefficient, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                coefficient (1,1) radia.ngsolve.CoefficientFunction
                options.Label (1,1) string = ...
                    "radia_matlab_coefficient_linear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.linear_form.create_from_coefficient', ...
                space.nativeHandle(), coefficient.nativeHandle(), ...
                char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.linear_form.info', nativeHandle);
            obj = radia.ngsolve.LinearForm(nativeHandle, info);
        end

        function obj = createBoundaryFromCoefficient(space, coefficient, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                coefficient (1,1) radia.ngsolve.CoefficientFunction
                options.Label (1,1) string = ...
                    "radia_matlab_boundary_coefficient_linear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.linear_form.create_boundary_from_coefficient', ...
                space.nativeHandle(), coefficient.nativeHandle(), ...
                char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.linear_form.info', nativeHandle);
            obj = radia.ngsolve.LinearForm(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex('ngsolve.linear_form.info', ...
                obj.NativeHandle);
        end

        function vector = vector(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.linear_form.vector', obj.NativeHandle);
            vector = radia.ngsolve.Vector.fromNativeHandle(nativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('ngsolve.linear_form.destroy', ...
                        obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:ngsolve:LinearFormDeleted", ...
                    "The native NGSolve LinearForm has been deleted.");
            end
        end
    end
end
