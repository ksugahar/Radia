classdef test_radiafield_mex < matlab.unittest.TestCase
    % Exercise the MATLAB bridge to the shared source-field cache and precision.
    properties (Access=private)
        Field
        MeshPath
    end
    methods (TestMethodSetup)
        function setupField(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'matlab')));
            radia.setup();
            radia.UtiDelAll();
            testCase.addTeardown(@() radia.UtiDelAll());
            radia.FldLenRndSw('off');
            testCase.addTeardown(@() radia.FldLenRndSw('on'));
            testCase.addTeardown(@() radia.FldCmpPrc('PrcB->0.0001,PrcA->0.001'));
            object = radia.ObjRecMag([-1,-1,-1], [0.1,0.1,0.1], [0,0,1]);
            testCase.Field = radia.RadiaField(object, 'b', Precision=2.5e-13);
            testCase.addTeardown(@() delete(testCase.Field));
            testCase.MeshPath = string(tempname('C:\temp')) + '.vol';
            lines = ["mesh3d"; "dimension"; "3"; "geomtype"; "0"; ...
                "facedescriptors"; "1"; "1 1 0 1 1"; "surfaceelements"; "4"; ...
                "1 1 1 0 3 1 2 3"; "1 1 1 0 3 1 4 2"; ...
                "1 1 1 0 3 2 4 3"; "1 1 1 0 3 3 4 1"; ...
                "volumeelements"; "1"; "1 4 1 2 3 4"; "points"; "4"; ...
                "0 0 0"; "1 0 0"; "0 1 0"; "0 0 1"; "pointelements"; "0"; ...
                "materials"; "1"; "1 air"; "bcnames"; "1"; "1 outer"; "endmesh"];
            writelines(lines, testCase.MeshPath);
            testCase.addTeardown(@() delete(testCase.MeshPath));
        end
    end
    methods (Test)
        function memoizationReusesValuesAndCanBeDisabled(testCase)
            field = testCase.Field;
            points = [0.15,0.20,0.10; 0.25,0.10,0.15];
            testCase.verifyFalse(field.fieldInfo().memoize);
            reference = field.evaluate(testCase.MeshPath, points);
            field.setMemoize(true);
            first = field.evaluate(testCase.MeshPath, points);
            second = field.evaluate(testCase.MeshPath, points);
            testCase.verifyTrue(field.fieldInfo().memoize);
            testCase.verifyEqual(second, first);
            testCase.verifyEqual(first, reference, RelTol=3e-13, AbsTol=2e-14);
            testCase.verifyEqual(field.cacheStats().size, 2);
            testCase.verifyGreaterThanOrEqual(field.cacheStats().hits, 2);
            field.setMemoize(false);
            testCase.verifyFalse(field.fieldInfo().memoize);
            field.clearCache();
            testCase.verifyEqual(field.cacheStats().size, 0);
        end
        function subMicroPrecisionIsPreserved(testCase)
            testCase.verifyEqual(testCase.Field.fieldInfo().precision, 2.5e-13);
            value = testCase.Field.evaluate(testCase.MeshPath, [0.15,0.20,0.10]);
            testCase.verifyTrue(all(isfinite(value), 'all'));
        end
    end
end
