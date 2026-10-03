function run_installed_wheel_acceptance(installedMatlabPath, evidencePath)
%RUN_INSTALLED_WHEEL_ACCEPTANCE Run the installed-wheel test and write its evidence.
%
% Every variable lives in this function workspace, so a borrowed MATLAB
% session keeps its base-workspace variables. The evidence file is closed on
% success and on failure.

arguments
    installedMatlabPath (1,1) string
    evidencePath (1,1) string
end

result = test_standalone_simulink(installedMatlabPath);
assert(result.ok, "radia:optuna:InstalledWheelAcceptance", ...
    "The installed-wheel Simulink test did not pass.");
fileId = fopen(evidencePath, "w");
assert(fileId >= 0, "radia:optuna:InstalledWheelEvidence", ...
    "Cannot write installed-wheel evidence: %s", evidencePath);
closeFile = onCleanup(@()fclose(fileId));
fprintf(fileId, "%s\n", jsonencode(result, PrettyPrint=true));
clear closeFile
end
