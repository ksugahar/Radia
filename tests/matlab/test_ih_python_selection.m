function tests = test_ih_python_selection
tests = functiontests(localfunctions);
end

function testSelectedInterpreterStateRoundTrip(testCase)
cleanup = onCleanup(@() radia.internal.selectedPythonExecutable("clear")); %#ok<NASGU>
chosen = "C:\temp\radia-selected-python.exe";
radia.internal.selectedPythonExecutable("set", chosen);
testCase.verifyEqual(radia.internal.selectedPythonExecutable("get"), chosen);
testCase.verifyEqual(radia.internal.resolveAssemblyPython(""), chosen);
override = "C:\temp\explicit-python.exe";
testCase.verifyEqual(radia.internal.resolveAssemblyPython(override), override);
end

function testSetupPublishesCachedInterpreter(testCase)
cleanup = onCleanup(@() radia.internal.selectedPythonExecutable("clear")); %#ok<NASGU>
python = string(getenv("RADIA_PYTHON_EXECUTABLE"));
testCase.assumeNotEmpty(python);
info = radia.setup(PythonExecutable=python, RequireMex=false, ...
    ConfigureSimulinkFileGeneration=false, Force=true);
testCase.verifyEqual( ...
    radia.internal.selectedPythonExecutable("get"), info.python_executable);
cached = radia.setup(PythonExecutable=python, RequireMex=false, ...
    ConfigureSimulinkFileGeneration=false);
testCase.verifyEqual( ...
    radia.internal.selectedPythonExecutable("get"), cached.python_executable);
end
