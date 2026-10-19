function value = pythonProcessPath(action)
%PYTHONPROCESSPATH Preserve the pre-MEX DLL path for external Python.

arguments
    action (1,1) string {mustBeMember(action, ["capture", "get"])}
end

persistent originalPath
if action == "capture" && isempty(originalPath)
    originalPath = string(getenv("PATH"));
end
if isempty(originalPath)
    value = string(getenv("PATH"));
else
    value = originalPath;
end

% A long-lived MATLAB session may have loaded an older Radia setup before
% this helper first captured PATH.  External CPython must never inherit the
% Netgen/NGSolve package DLL directories even in that lifecycle ordering.
if ispc
    entries = split(value, pathsep);
    normalized = lower(replace(entries, "/", "\"));
    unsafe = contains(normalized, "\lib\site-packages\ngsolve_openblas") | ...
        contains(normalized, "\lib\site-packages\netgen");
    value = strjoin(entries(~unsafe & entries ~= ""), pathsep);
end
end
