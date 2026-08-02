function tests = test_python_fallback
% Verify the explicit in-process Python fallback contract and conversion.
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
end

function testAcousticFallbackReportsBackendAndDll(testCase)
zeta = [0.1+0.2i, -0.4+0.3i];
result = radia.python.acoustic("cq", "bdf_delta", {zeta}, ...
    Keywords=struct("method", "BDF2"));
expected = 1.5 - 2*zeta + 0.5*zeta.^2;
verifyEqual(testCase, result.value, expected, RelTol=0, AbsTol=2e-15);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.python.execution_mode, "InProcess");
verifyTrue(testCase, endsWith(lower(result.python.library), "python312.dll"));
end

function testFamilyWrapperRejectsWrongModule(testCase)
call = @() radia.python.acoustic("not_acoustic", "run", {});
verifyError(testCase, call, "MATLAB:validators:mustBeMember");
end

function testNumpyScalarConvertsToMatlabNumeric(testCase)
result = radia.internal.callPython("numpy", "float64", {1.25});
verifyEqual(testCase, result.value, 1.25, "AbsTol", 0);
verifyClass(testCase, result.value, "double");
end

function testIsochronousTopoptHasNamedFallback(testCase)
result = radia.python.isochronousTopopt( ...
    "density_to_s", {[0, 1], 100});
verifyEqual(testCase, result.value, [1e6, 0.01], "RelTol", 1e-14);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.python.execution_mode, "InProcess");
end

function testHarmonicBalanceHasNamedFallback(testCase)
result = radia.python.harmonicBalance("periodic_phase", {8});
verifyEqual(testCase, result.value, (0:7) * pi / 4, "AbsTol", 2e-15);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.harmonic_balance");
end

function testXsuiteBridgeHasNamedFallback(testCase)
boundary = radia.python.xsuiteBridge( ...
    "AxisAlignedBox", {[-1, -1, -1], [1, 1, 1]});
x = [0; 0; 2];
y = zeros(3, 1);
z = zeros(3, 1);
result = radia.python.xsuiteBridge( ...
    "first_box_exit_events", {x, y, z, boundary.value});
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.xsuite_bridge");
verifyEqual(testCase, result.value.particle_index, 0);
verifyEqual(testCase, result.value.step_index, 2);
verifyEqual(testCase, result.value.position_m, {2, 0, 0});
end

function testExternalPythonUsesChildSafeDllPathAndRestoresMatlabPath(testCase)
if ~ispc
    return
end
radia.setup(Force=true);
activePath = string(getenv("PATH"));
[status, childPath] = radia.internal.runPythonProcess("echo %PATH%");
verifyEqual(testCase, status, 0);
normalized = lower(replace(string(strtrim(childPath)), "/", "\"));
verifyFalse(testCase, contains(normalized, ...
    "\lib\site-packages\ngsolve_openblas"));
verifyFalse(testCase, contains(normalized, ...
    "\lib\site-packages\netgen"));
verifyEqual(testCase, string(getenv("PATH")), activePath);
end
