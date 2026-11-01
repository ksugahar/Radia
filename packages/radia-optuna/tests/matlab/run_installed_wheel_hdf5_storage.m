function run_installed_wheel_hdf5_storage(installedMatlabPath, testFile, evidencePath)
%RUN_INSTALLED_WHEEL_HDF5_STORAGE Run the HDF5 storage tests on an installed wheel.
%
% The caller puts installedMatlabPath first on the path. Every entry point
% the tests use must resolve inside the installed wheel, and the evidence
% records the resolved files, the MEX build identity and each test result.
% The evidence is written before the pass/fail assertion, so a failed run
% still leaves its receipt.

arguments
    installedMatlabPath (1,1) string
    testFile (1,1) string
    evidencePath (1,1) string
end

names = ["radia.optuna.Study", "radia.optuna.OptimizationSession", ...
    "radia.optuna.internal.TypedHdf5", "optuna_mex"];
resolved = strings(size(names));
for k = 1:numel(names)
    resolved(k) = string(which(names(k)));
    assert(startsWith(lower(resolved(k)), lower(installedMatlabPath)), ...
        "radia:optuna:InstalledWheelPath", ...
        "%s resolved outside the installed wheel: %s", names(k), resolved(k));
end
info = optuna_mex('api.info');

previous = getenv("RADIA_OPTUNA_MATLAB_PATH");
restore = onCleanup(@() setenv("RADIA_OPTUNA_MATLAB_PATH", previous));
setenv("RADIA_OPTUNA_MATLAB_PATH", installedMatlabPath);
results = runtests(char(testFile));
clear restore

receipt = struct();
receipt.schema = "radia.optuna.installed-wheel-hdf5-storage.v1";
receipt.matlab_release = string(version("-release"));
receipt.installed_matlab_path = installedMatlabPath;
receipt.test_file = testFile;
receipt.resolved = cell2struct(cellstr(resolved(:)), ...
    matlab.lang.makeValidName(cellstr(names(:))), 1);
receipt.mex_info = info;
receipt.tests = struct("name", {results.Name}, ...
    "passed", {results.Passed}, "failed", {results.Failed}, ...
    "incomplete", {results.Incomplete}, "duration_s", {results.Duration});
receipt.ok = ~isempty(results) && all([results.Passed]);

fileId = fopen(evidencePath, "w");
assert(fileId >= 0, "radia:optuna:InstalledWheelEvidence", ...
    "Cannot write installed-wheel evidence: %s", evidencePath);
closeFile = onCleanup(@()fclose(fileId));
fprintf(fileId, "%s\n", jsonencode(receipt, PrettyPrint=true));
clear closeFile
assert(receipt.ok, "radia:optuna:InstalledWheelHdf5Storage", ...
    "%d of %d installed-wheel HDF5 storage tests passed.", ...
    nnz([results.Passed]), numel(results));
end
