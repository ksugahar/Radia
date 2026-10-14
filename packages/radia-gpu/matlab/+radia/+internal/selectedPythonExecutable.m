function executable = selectedPythonExecutable(action, value)
%SELECTEDPYTHONEXECUTABLE Process-local Python selected by radia.setup.
%   This state lets explicit assembly boundaries reuse the interpreter whose
%   NGSolve ABI radia.setup already validated, without searching PATH again.

arguments
    action (1,1) string {mustBeMember(action, ["get", "set", "clear"])}
    value (1,1) string = ""
end
persistent selected
if action == "set"
    if strlength(strtrim(value)) == 0
        error("radia:setup:Python", "Selected Python executable is empty.");
    end
    selected = value;
elseif action == "clear"
    selected = "";
end
if isempty(selected)
    executable = "";
else
    executable = selected;
end
end
