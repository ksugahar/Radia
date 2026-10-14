classdef test_ih_python < matlab.unittest.TestCase
    % Public IH batch adapters: keyword arguments, records and Python objects.
    methods (TestClassSetup)
        function addMatlabPackage(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'matlab')));
        end
    end
    methods (Test)
        function thermalPowerTransfer(testCase)
            result = radia.python.ihThermal("power_balance", {101, 100, 0.02}, ...
                Keywords=struct("what", "MATLAB transfer check"));
            testCase.verifyEqual(result.backend, "python-fallback");
            testCase.verifyEqual(result.value.relative_error, 0.01, AbsTol=1e-15);
            testCase.verifyEqual(result.value.target_W, 101);
        end
        function thermalExposureConstraint(testCase)
            exposure = struct("T_max_C", 850, "thresholds", ...
                {{struct("T_C", 800, "volume_m3", 2e-6)}});
            result = radia.python.ihThermalPost("limit_check", {exposure, 800});
            testCase.verifyTrue(result.value.exceeded);
            testCase.verifyEqual(result.value.excess_C, 50);
            testCase.verifyEqual(result.value.volume_above_m3, 2e-6);
        end
        function coupledMaterialTable(testCase)
            result = radia.python.ihAxisymCoupled("EMMaterialTable", ...
                {[20, 100], [6e7, 5e7], [100, 10]}, Keywords=struct("source", "MATLAB test"));
            values = cell(result.value.evaluate(60));
            testCase.verifyEqual(double(values{1}), 5.5e7, RelTol=1e-14);
            testCase.verifyEqual(double(values{2}), 55, AbsTol=1e-13);
            testCase.verifyEqual(string(result.value.source), "MATLAB test");
        end
    end
end
