classdef PlanarFieldEvaluator < handle
    %PLANARFIELDEVALUATOR Persistent planar charge-cloud field evaluator.

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = PlanarFieldEvaluator(handle)
            obj.NativeHandle = uint64(handle);
        end
    end

    methods (Static)
        function obj = create(positions, strengths, imageMasks, imageSigns)
            arguments
                positions double
                strengths double
                imageMasks = int32.empty
                imageSigns double = zeros(0,1)
            end
            handle = radia.internal.callMex('hdiv.planar_evaluator.create', ...
                positions, strengths, int32(imageMasks), double(imageSigns));
            obj = radia.PlanarFieldEvaluator(handle);
        end

        function obj = fromChargeGram(chargeGram, magnetization)
            arguments
                chargeGram (1,1) radia.HACApKChargeGram
                magnetization double
            end
            handle = radia.internal.callMex( ...
                'hacapk.charge_gram.create_planar_field_evaluator', ...
                chargeGram.nativeHandle(), magnetization);
            obj = radia.PlanarFieldEvaluator(handle);
        end
    end

    methods
        function value = field(obj, points)
            obj.assertAlive();
            value = radia.internal.callMex('hdiv.planar_evaluator.field', ...
                obj.NativeHandle, points);
        end

        function value = az(obj, points)
            obj.assertAlive();
            value = radia.internal.callMex('hdiv.planar_evaluator.az', ...
                obj.NativeHandle, points);
        end

        function value = stats(obj)
            obj.assertAlive();
            value = radia.internal.callMex('hdiv.planar_evaluator.stats', obj.NativeHandle);
        end

        function coefficient = coefficientFunction(obj, options)
            arguments
                obj (1,1) radia.PlanarFieldEvaluator
                options.SourceAngle (1,1) double = 0
                options.TargetAngle (1,1) double = 0
                options.Center (1,2) double = [0,0]
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'hdiv.planar_evaluator.as_coefficient', obj.NativeHandle, ...
                options.SourceAngle, options.TargetAngle, ...
                options.Center(1), options.Center(2));
            coefficient = radia.ngsolve.CoefficientFunction.fromNativeHandle( ...
                nativeHandle);
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('hdiv.planar_evaluator.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:PlanarFieldEvaluator:Deleted", ...
                    "The native planar field evaluator has been deleted.");
            end
        end
    end
end
