function tests = test_picard_acceleration_fallback
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root, "matlab"));
end

function testStatefulAcceleratorMatchesPython(testCase)
result = radia.python.picardAcceleration("ConstrainedAndersonAccelerator", ...
    {2, 0.3}, Keywords=struct("lower", 1, "upper", 100));
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.picard_acceleration");
actual = result.value;
module = py.importlib.import_module("radia.picard_acceleration");
oracle = module.ConstrainedAndersonAccelerator(int32(2), 0.3, ...
    pyargs("lower", 1, "upper", 100));
x = py.numpy.asarray([2, 8]);
for k = 1:4
    target = py.numpy.asarray([3+k/10, 5+k/10]);
    next = actual.step(x, target);
    expected = oracle.step(x, target);
    verifyEqual(testCase, double(next), double(expected), AbsTol=0);
    x = next;
end
actual.reset();
oracle.reset();
verifyEqual(testCase, double(actual.step(x, target)), ...
    double(oracle.step(x, target)), AbsTol=0);
end

function testContractionDiagnosticMatchesPython(testCase)
values = [1, .5, .25, .125];
result = radia.python.picardAcceleration("estimate_contraction_rate", {values});
module = py.importlib.import_module("radia.picard_acceleration");
expected = module.estimate_contraction_rate(py.numpy.asarray(values));
verifyEqual(testCase, result.value, double(expected), AbsTol=0);
end
