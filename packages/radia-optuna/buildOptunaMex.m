function binary = buildOptunaMex(python,buildDirectory,pythonVerification)
%BUILDOPTUNAMEX Link a plain MSVC MEX to the Python-verified nvcc archive.
% Build and verify the package Python kernel target first. No compiler fallback.
arguments
    python (1,1) string
    buildDirectory (1,1) string
    pythonVerification (1,1) string
end
assert(ispc && isfile(python) && isfile(pythonVerification), ...
    'radia:optuna:Build','Explicit Python and verification receipt are required.');
package=fileparts(mfilename('fullpath'));
script=fullfile(package,'tools','build_native.ps1');
command=sprintf('pwsh -NoProfile -File "%s" -Stage matlab -Python "%s" -MatlabRoot "%s" -BuildDirectory "%s" -PythonVerification "%s"', ...
    script,python,matlabroot,buildDirectory,pythonVerification);
[status,output]=system(command);fprintf('%s',output);
assert(status==0,'radia:optuna:Build','Plain MEX/shared-library build failed.');
binary=fullfile(package,'matlab',"optuna_mex."+mexext);
assert(isfile(binary),'radia:optuna:Build','Plain MEX target is missing.');
end
