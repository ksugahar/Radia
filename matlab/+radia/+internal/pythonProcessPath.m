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
end
