function tests = test_installed_wheel_hdf5_gate
%TEST_INSTALLED_WHEEL_HDF5_GATE Negative checks for the installed-wheel gate.
%
% Guards: the inventory accepts only the exact expected named set with
% every test passed, and path containment rejects same-prefix siblings
% while accepting either separator convention.
% Assumes: no MEX, GPU or wheel; pure gate logic.
tests = functiontests(localfunctions);
end

function testExactInventoryPasses(testCase)
[ok, report] = installed_wheel_inventory(fakeResults(expected()), expected());
verifyTrue(testCase, ok, jsonencode(report));
end

function testOmittedTestsAreRejected(testCase)
names = expected();
[ok, report] = installed_wheel_inventory(fakeResults(names(1)), names);
verifyFalse(testCase, ok);
verifyEqual(testCase, report.executed_count, 1);
verifyEqual(testCase, numel(report.missing), 7);
end

function testExtraAndDuplicatedTestsAreRejected(testCase)
names = expected();
[ok, report] = installed_wheel_inventory(fakeResults([names, "x/extra"]), names);
verifyFalse(testCase, ok);
verifyEqual(testCase, report.extra, "x/extra");
[ok, report] = installed_wheel_inventory(fakeResults([names, names(2)]), names);
verifyFalse(testCase, ok);
verifyEqual(testCase, report.duplicated, names(2));
end

function testFailedAndIncompleteTestsAreRejected(testCase)
names = expected();
results = fakeResults(names);
results(3).Passed = false;
results(3).Failed = true;
[ok, report] = installed_wheel_inventory(results, names);
verifyFalse(testCase, ok);
verifyEqual(testCase, report.not_passed, names(3));
results = fakeResults(names);
results(5).Passed = false;
results(5).Incomplete = true;
[ok, report] = installed_wheel_inventory(results, names);
verifyFalse(testCase, ok);
verifyEqual(testCase, report.not_passed, names(5));
end

function testPathContainment(testCase)
root = "C:/temp/run03/layout";
verifyTrue(testCase, installed_wheel_path_inside( ...
    "C:\temp\run03\layout\+radia\+optuna\Study.m", root));
verifyTrue(testCase, installed_wheel_path_inside( ...
    "C:/temp/run03/layout/optuna_mex.mexw64", root + "/"));
verifyFalse(testCase, installed_wheel_path_inside( ...
    "C:/temp/run03/layout-sibling/+radia/+optuna/Study.m", root));
verifyFalse(testCase, installed_wheel_path_inside(root, root));
verifyFalse(testCase, installed_wheel_path_inside("", root));
verifyError(testCase, @() installed_wheel_path_inside( ...
    "C:/temp/run03/layout/../other/Study.m", root), ...
    "radia:optuna:InstalledWheelPath");
end

function names = expected()
names = "test_optuna_hdf5_storage/t" + string(1:8);
end

function results = fakeResults(names)
results = struct("Name", cellstr(names), "Passed", true, ...
    "Failed", false, "Incomplete", false);
end
