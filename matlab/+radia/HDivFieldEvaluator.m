classdef HDivFieldEvaluator < handle
    %HDIVFIELDEVALUATOR Persistent direct/tree HDiv field evaluator.

    properties (Access=private)
        NativeHandle uint64 = uint64(0)
    end

    methods (Access=private)
        function obj = HDivFieldEvaluator(handle)
            obj.NativeHandle = uint64(handle);
        end
    end

    methods (Static)
        function obj = fromTet(volume, surface, imageMasks, imageSigns, options)
            arguments
                volume double
                surface double
                imageMasks = int32.empty
                imageSigns double = zeros(0,1)
                options.LeafSize (1,1) double = 32
                options.Theta (1,1) double = 0.05
                options.TreeMinSources (1,1) double = 256
                options.AutoMinWork (1,1) double = 500000000
                options.TreeRelativeTolerance (1,1) double = 1e-5
                options.ProbeCount (1,1) double = 16
            end
            handle = radia.internal.callMex('hdiv.field_evaluator.from_tet', ...
                volume, surface, int32(imageMasks), double(imageSigns), ...
                options.LeafSize, options.Theta, options.TreeMinSources, ...
                options.AutoMinWork, options.TreeRelativeTolerance, options.ProbeCount);
            obj = radia.HDivFieldEvaluator(handle);
        end

        function obj = fromCloud(xyz, strength, imageMasks, imageSigns, options)
            arguments
                xyz double
                strength double
                imageMasks = int32.empty
                imageSigns double = zeros(0,1)
                options.LeafSize (1,1) double = 32
                options.Theta (1,1) double = 0.05
                options.TreeMinSources (1,1) double = 256
                options.AutoMinWork (1,1) double = 500000000
                options.TreeRelativeTolerance (1,1) double = 1e-5
                options.ProbeCount (1,1) double = 16
            end
            handle = radia.internal.callMex('hdiv.field_evaluator.from_cloud', ...
                xyz, strength, int32(imageMasks), double(imageSigns), ...
                options.LeafSize, options.Theta, options.TreeMinSources, ...
                options.AutoMinWork, options.TreeRelativeTolerance, options.ProbeCount);
            obj = radia.HDivFieldEvaluator(handle);
        end

        function obj = fromCurvedTet(volume, surface, gaussPoints, gaussWeights, ...
                imageMasks, imageSigns, options)
            arguments
                volume double
                surface double
                gaussPoints double
                gaussWeights double
                imageMasks = int32.empty
                imageSigns double = zeros(0,1)
                options.LeafSize (1,1) double = 32
                options.Theta (1,1) double = 0.05
                options.TreeMinSources (1,1) double = 256
                options.AutoMinWork (1,1) double = 500000000
                options.TreeRelativeTolerance (1,1) double = 1e-5
                options.ProbeCount (1,1) double = 16
            end
            handle = radia.internal.callMex('hdiv.field_evaluator.from_curved_tet', ...
                volume, surface, gaussPoints, gaussWeights, int32(imageMasks), ...
                double(imageSigns), options.LeafSize, options.Theta, ...
                options.TreeMinSources, options.AutoMinWork, ...
                options.TreeRelativeTolerance, options.ProbeCount);
            obj = radia.HDivFieldEvaluator(handle);
        end

        function obj = fromChargeGram(chargeGram, magnetization, options)
            arguments
                chargeGram (1,1) radia.HACApKChargeGram
                magnetization double
                options.LeafSize (1,1) double = 32
                options.Theta (1,1) double = 0.05
                options.TreeMinSources (1,1) double = 256
                options.AutoMinWork (1,1) double = 500000000
                options.TreeRelativeTolerance (1,1) double = 1e-5
                options.ProbeCount (1,1) double = 16
            end
            handle = radia.internal.callMex('hacapk.charge_gram.create_field_evaluator', ...
                chargeGram.nativeHandle(), magnetization, options.LeafSize, options.Theta, ...
                options.TreeMinSources, options.AutoMinWork, ...
                options.TreeRelativeTolerance, options.ProbeCount);
            obj = radia.HDivFieldEvaluator(handle);
        end
    end

    methods
        function value = field(obj, observations, options)
            arguments
                obj (1,1) radia.HDivFieldEvaluator
                observations double
                options.Algorithm (1,1) string = "auto"
            end
            obj.assertAlive();
            value = radia.internal.callMex('hdiv.field_evaluator.field', ...
                obj.NativeHandle, observations, char(options.Algorithm));
        end

        function value = candidateAlgorithm(obj, nObservations)
            obj.assertAlive();
            value = string(radia.internal.callMex( ...
                'hdiv.field_evaluator.candidate_algorithm', obj.NativeHandle, nObservations));
        end

        function value = lastAlgorithm(obj)
            obj.assertAlive();
            value = string(radia.internal.callMex( ...
                'hdiv.field_evaluator.last_algorithm', obj.NativeHandle));
        end

        function value = stats(obj)
            obj.assertAlive();
            value = radia.internal.callMex('hdiv.field_evaluator.stats', obj.NativeHandle);
        end

        function coefficient = coefficientFunction(obj, options)
            arguments
                obj (1,1) radia.HDivFieldEvaluator
                options.Algorithm (1,1) string = "direct"
            end
            obj.assertAlive();
            nativeHandle = radia.internal.callMex( ...
                'hdiv.field_evaluator.as_coefficient', obj.NativeHandle, ...
                char(options.Algorithm));
            coefficient = radia.ngsolve.CoefficientFunction.fromNativeHandle( ...
                nativeHandle);
        end

        function handle = nativeHandle(obj)
            obj.assertAlive();
            handle = obj.NativeHandle;
        end

        function delete(obj)
            if obj.NativeHandle ~= 0
                try
                    radia.internal.callMex('hdiv.field_evaluator.destroy', obj.NativeHandle);
                catch
                end
                obj.NativeHandle = uint64(0);
            end
        end
    end

    methods (Access=private)
        function assertAlive(obj)
            if obj.NativeHandle == 0
                error("radia:HDivFieldEvaluator:Deleted", ...
                    "The native HDiv field evaluator has been deleted.");
            end
        end
    end
end
