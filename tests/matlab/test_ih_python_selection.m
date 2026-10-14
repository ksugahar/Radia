classdef test_ih_python_selection < matlab.unittest.TestCase
    %TEST_IH_PYTHON_SELECTION Interpreter selection and setup cache contract.

    methods (TestMethodSetup)
        function preserveSelectedInterpreter(testCase)
            previous = radia.internal.selectedPythonExecutable("get");
            testCase.addTeardown(@() restoreSelectedInterpreter(previous));
        end
    end

    methods (Test)
        function selectedInterpreterStateRoundTrip(testCase)
            chosen = "C:\temp\radia-selected-python.exe";
            override = "C:\temp\explicit-python.exe";

            radia.internal.selectedPythonExecutable("set", chosen);
            selected = radia.internal.selectedPythonExecutable("get");
            resolved = radia.internal.resolveAssemblyPython("");
            explicit = radia.internal.resolveAssemblyPython(override);

            testCase.verifyEqual(selected, chosen);
            testCase.verifyEqual(resolved, chosen);
            testCase.verifyEqual(explicit, override);
        end

        function setupPublishesCachedInterpreter(testCase)
            python = string(getenv("RADIA_PYTHON_EXECUTABLE"));
            testCase.assumeNotEmpty(python);

            info = radia.setup(PythonExecutable=python, RequireMex=false, ...
                ConfigureSimulinkFileGeneration=false, Force=true);
            selected = radia.internal.selectedPythonExecutable("get");
            cached = radia.setup(PythonExecutable=python, RequireMex=false, ...
                ConfigureSimulinkFileGeneration=false);
            selectedCached = radia.internal.selectedPythonExecutable("get");

            testCase.verifyEqual(selected, info.python_executable);
            testCase.verifyEqual(selectedCached, cached.python_executable);
        end
    end
end

function restoreSelectedInterpreter(previous)
if strlength(previous) == 0
    radia.internal.selectedPythonExecutable("clear");
else
    radia.internal.selectedPythonExecutable("set", previous);
end
end
