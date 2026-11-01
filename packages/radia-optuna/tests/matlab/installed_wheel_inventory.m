function [ok, report] = installed_wheel_inventory(results, expected)
%INSTALLED_WHEEL_INVENTORY Check that exactly the expected tests ran and passed.
%
% results is a runtests result array (or a struct array with Name, Passed,
% Failed and Incomplete). ok is true only when the executed names equal the
% expected set without duplicates and every test passed, did not fail and
% was not incomplete. report lists what is missing, extra, duplicated or
% did not pass.

arguments
    results
    expected (1,:) string
end
names = string({results.Name});
passed = logical([results.Passed]);
failed = logical([results.Failed]);
incomplete = logical([results.Incomplete]);
[unique_names, first] = unique(names);
report = struct();
report.expected_count = numel(expected);
report.executed_count = numel(names);
report.missing = setdiff(expected, names);
report.extra = setdiff(names, expected);
report.duplicated = unique(names(setdiff(1:numel(names), first)));
report.not_passed = names(~passed | failed | incomplete);
ok = numel(unique(expected)) == numel(expected) && ...
    numel(unique_names) == numel(names) && ...
    isempty(report.missing) && isempty(report.extra) && ...
    isempty(report.not_passed) && numel(names) == numel(expected);
end
