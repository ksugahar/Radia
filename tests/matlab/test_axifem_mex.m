function tests = test_axifem_mex
% Cross-language parity tests for the shared-native axifem Q1/Q2 boundary.
tests = functiontests(localfunctions);
end

function testQ2InteriorAndAxisElementsMatchPybind(testCase)
reference = testCase.TestData.reference;
verifyQ2Case(testCase, reference.q2_interior_input, ...
    reference.q2_interior_stiffness, reference.q2_interior_sigma_mass);
actual = verifyQ2Case(testCase, reference.q2_axis_input, ...
    reference.q2_axis_stiffness, reference.q2_axis_sigma_mass);
verifyTrue(testCase, actual.axis_touching);
axisDofs = [1, 4, 8];
verifyEqual(testCase, actual.stiffness(axisDofs, :), zeros(3, 9), AbsTol=0);
verifyEqual(testCase, actual.stiffness(:, axisDofs), zeros(9, 3), AbsTol=0);
verifyEqual(testCase, actual.sigma_mass(axisDofs, :), zeros(3, 9), AbsTol=0);
verifyEqual(testCase, actual.sigma_mass(:, axisDofs), zeros(9, 3), AbsTol=0);
end

function setupOnce(testCase)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
radia.setup(Force=true);
[path, reference] = loadPythonReference(repoRoot);
testCase.TestData.referencePath = path;
testCase.TestData.reference = reference;
end

function teardownOnce(testCase)
if isfile(testCase.TestData.referencePath)
    delete(testCase.TestData.referencePath);
end
end

function testInteriorAndAxisElementsMatchPybind(testCase)
reference = testCase.TestData.reference;
verifyCase(testCase, reference.interior_input, ...
    reference.interior_stiffness, reference.interior_sigma_mass);
verifyCase(testCase, reference.axis_input, ...
    reference.axis_stiffness, reference.axis_sigma_mass);
verifyCase(testCase, reference.thin_input, ...
    reference.thin_stiffness, reference.thin_sigma_mass);
verifyCase(testCase, reference.insulator_input, ...
    reference.insulator_stiffness, reference.insulator_sigma_mass);
verifyEqual(testCase, reference.insulator_sigma_mass, zeros(4), AbsTol=0);
end

function testCommandAndMetadataAreExplicit(testCase)
commands = string(radia.internal.callMex('api.commands'));
verifyTrue(testCase, any(commands == "axifem.q1_magnetic_element_matrices"));
verifyTrue(testCase, any(commands == "axifem.q2_magnetic_element_matrices"));
values = testCase.TestData.reference.interior_input;
result = radia.axifem.q1MagneticElementMatrices(values(1), values(2), ...
    values(3), values(4), values(5), values(6));
verifyEqual(testCase, string(result.backend), "native-mex");
verifyEqual(testCase, string(result.dof_convention), "nodal A_phi (V-DOF)");
verifyEqual(testCase, string(result.node_order), ...
    "(ra,za),(rb,za),(rb,zb),(ra,zb)");

q2 = radia.axifem.q2MagneticElementMatrices(values(1), values(2), ...
    values(3), values(4), values(5), values(6));
verifyEqual(testCase, string(q2.backend), "native-mex");
verifyEqual(testCase, string(q2.dof_convention), "nodal A_phi (V-DOF)");
verifyEqual(testCase, size(q2.stiffness), [9, 9]);
end

function actual = verifyQ2Case(testCase, values, expectedK, expectedM)
actual = radia.axifem.q2MagneticElementMatrices(values(1), values(2), ...
    values(3), values(4), values(5), values(6));
verifyEqual(testCase, actual.stiffness, expectedK, ...
    RelTol=2e-12, AbsTol=2e-10*max(1, max(abs(expectedK), [], "all")));
verifyEqual(testCase, actual.sigma_mass, expectedM, ...
    RelTol=2e-12, AbsTol=2e-10*max(1, max(abs(expectedM), [], "all")));
end

function testInvalidGeometryFailsLoudly(testCase)
call = @() radia_mex('axifem.q1_magnetic_element_matrices', ...
    2.0e-3, 1.0e-3, 0.0, 1.0e-3, 4*pi*1e-7, 5.8e7);
verifyError(testCase, call, "radia:mex:Exception");
end

function verifyCase(testCase, values, expectedK, expectedM)
actual = radia.axifem.q1MagneticElementMatrices(values(1), values(2), ...
    values(3), values(4), values(5), values(6));
verifyEqual(testCase, actual.stiffness, expectedK, ...
    RelTol=2e-13, AbsTol=2e-12*max(1, max(abs(expectedK), [], "all")));
verifyEqual(testCase, actual.sigma_mass, expectedM, ...
    RelTol=2e-13, AbsTol=2e-12*max(1, max(abs(expectedM), [], "all")));
end

function [path, reference] = loadPythonReference(repoRoot)
path = string(tempname("C:\temp")) + ".mat";
pythonExecutable = string(getenv("RADIA_PYTHON_EXECUTABLE"));
if strlength(pythonExecutable) == 0
    pythonExecutable = "python";
end
helper = fullfile(repoRoot, "tests", "matlab", "axifem_python_reference.py");
command = quoteCommandArgument(pythonExecutable) + " " + ...
    quoteCommandArgument(helper) + " " + quoteCommandArgument(path);
[status, output] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:test:PythonReference", ...
        "Python axifem reference failed (%d): %s", status, output);
end
reference = load(path);
end

function value = quoteCommandArgument(value)
value = '"' + replace(string(value), '"', '""') + '"';
end
