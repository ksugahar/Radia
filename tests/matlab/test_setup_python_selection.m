classdef test_setup_python_selection < matlab.unittest.TestCase
    properties
        Python
    end
    methods (TestMethodSetup)
        function configure(testCase)
            testCase.Python = string(getenv('RADIA_PYTHON_EXECUTABLE'));
            testCase.assertNotEmpty(char(testCase.Python), 'Provide the accepted Python runtime');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'RADIA_PYTHON_EXECUTABLE', 'C:/temp/radia-nonexistent-test-python.exe'));
            testCase.assertFalse(isfile(getenv('RADIA_PYTHON_EXECUTABLE')));
        end
    end
    methods (Test)
        function environmentSelectsDefault(testCase)
            testCase.verifyError(@() radia.setup(Force=true, ...
                ConfigureSimulinkFileGeneration=false, RequireMex=false), 'radia:setup:Python');
        end
        function explicitChoiceOverridesEnvironment(testCase)
            info = radia.setup(PythonExecutable=testCase.Python, Force=true, ...
                ConfigureSimulinkFileGeneration=false, RequireMex=false);
            testCase.verifyEqual(info.python_executable, testCase.Python);
        end
    end
end

