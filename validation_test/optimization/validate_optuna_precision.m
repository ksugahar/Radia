function report=validate_optuna_precision(outputPath)
%VALIDATE_OPTUNA_PRECISION Run binary64 diagnostics with upstream test gates.
% Configure pyenv to the pinned oracle Python before calling this function.
arguments
    outputPath (1,1) string
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
previous=getenv('RADIA_OPTUNA_PRECISION_OUTPUT');
cleanup=onCleanup(@()setenv('RADIA_OPTUNA_PRECISION_OUTPUT',previous));
setenv('RADIA_OPTUNA_PRECISION_OUTPUT',outputPath);
diagnostics=runtests(fullfile(root,'tests','matlab','test_optuna_table.m'), ...
    Name='*FloatPrecisionDiagnosticContract*');
assert(numel(diagnostics)==1 && all([diagnostics.Passed]), ...
    'Precision diagnostic self-test missing or failed.');
results=runtests(fullfile(root,'tests','matlab','test_optuna_upstream_oracle.m'));
report=jsondecode(fileread(outputPath));
report.tests_passed=sum([results.Passed]);
report.tests_failed=sum([results.Failed]);
report.tests_incomplete=sum([results.Incomplete]);
report.diagnostic_self_test_passed=all([diagnostics.Passed]);
fid=fopen(outputPath,'w'); assert(fid>=0);
fileCleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
assert(~isempty(results) && all([results.Passed]), ...
    'Upstream precision comparison missing or failed.');
end
