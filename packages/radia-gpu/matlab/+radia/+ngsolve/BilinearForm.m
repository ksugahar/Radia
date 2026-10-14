classdef BilinearForm < handle
    %BILINEARFORM Persistent assembled native NGSolve BilinearForm.
    %   Built-in real integrators can use scalar CoefficientFunction weights.

    properties (SetAccess=private)
        Space string = ""
        Form string = ""
        Label string = ""
    end

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = BilinearForm(nativeHandle, info)
            obj.NativeHandle = uint64(nativeHandle);
            obj.Space = string(info.space);
            obj.Form = string(info.form);
            obj.Label = string(info.label);
        end
    end

    methods (Static)
        function obj = create(space, form, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                form (1,1) string
                options.Label (1,1) string = "radia_matlab_bilinear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.bilinear_form.create', space.nativeHandle(), ...
                char(form), char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.bilinear_form.info', nativeHandle);
            obj = radia.ngsolve.BilinearForm(nativeHandle, info);
        end

        function obj = createFromCoefficient(space, form, coefficient, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                form (1,1) string
                coefficient (1,1) radia.ngsolve.CoefficientFunction
                options.Label (1,1) string = ...
                    "radia_matlab_coefficient_bilinear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.bilinear_form.create_from_coefficient', ...
                space.nativeHandle(), char(form), coefficient.nativeHandle(), ...
                char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.bilinear_form.info', nativeHandle);
            obj = radia.ngsolve.BilinearForm(nativeHandle, info);
        end

        function obj = createBoundaryFromCoefficient(space, coefficient, options)
            arguments
                space (1,1) radia.ngsolve.FESpace
                coefficient (1,1) radia.ngsolve.CoefficientFunction
                options.Label (1,1) string = ...
                    "radia_matlab_boundary_coefficient_bilinear_form"
            end
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.bilinear_form.create_boundary_from_coefficient', ...
                space.nativeHandle(), coefficient.nativeHandle(), ...
                char(options.Label));
            info = radia.internal.callMex( ...
                'ngsolve.bilinear_form.info', nativeHandle);
            obj = radia.ngsolve.BilinearForm(nativeHandle, info);
        end
    end

    methods
        function value = info(obj)
            obj.assertAlive();
            value = radia.internal.callMex( ...
                'ngsolve.bilinear_form.info', obj.NativeHandle);
        end

        function matrix = matrix(obj)
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'ngsolve.bilinear_form.matrix', obj.NativeHandle);
            matrix = radia.ngsolve.Matrix.fromNativeHandle(nativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex( ...
                        'ngsolve.bilinear_form.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:ngsolve:BilinearFormDeleted", ...
                    "The native NGSolve BilinearForm has been deleted.");
            end
        end
    end
end
