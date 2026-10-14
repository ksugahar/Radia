function [status, output] = runPythonProcess(command)
%RUNPYTHONPROCESS Run an external Python command with a child-safe DLL path.
% radia.setup prepends Netgen and NGSolve directories so MATLAB can resolve
% radia_mex dependencies. A fresh Python process configures those directories
% itself; inheriting the MATLAB-only entries can load a conflicting BLAS DLL.

arguments
    command (1,1) string
end

radia.internal.pythonProcessPath("capture");
activePath = string(getenv("PATH"));
if ~ispc
    [status, output] = system(command);
    return
end

childPath = radia.internal.pythonProcessPath("get");
cleanup = onCleanup(@() setenv("PATH", char(activePath)));
setenv("PATH", char(childPath));
[status, output] = system(command);
clear cleanup
end
