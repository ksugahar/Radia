function executable = findExecutable(options)
%FINDEXECUTABLE Locate a user-installed LTspice executable.
arguments
    options.Executable (1,1) string = ""
end

if strlength(options.Executable) > 0
    if ~isfile(options.Executable)
        error("radia:ltspice:ExecutableNotFound", ...
            "LTspice executable does not exist: %s", options.Executable);
    end
    executable = options.Executable;
    return
end

candidates = strings(0, 1);
if ispc
    candidates = [
        "C:\Program Files\ADI\LTspice\LTspice.exe"
        string(getenv("LOCALAPPDATA")) + "\Programs\ADI\LTspice\LTspice.exe"
    ];
elseif ismac
    candidates = [
        "/Applications/LTspice.app/Contents/MacOS/LTspice"
        string(getenv("HOME")) + "/Applications/LTspice.app/Contents/MacOS/LTspice"
    ];
end

match = candidates(isfile(candidates));
if isempty(match)
    error("radia:ltspice:ExecutableNotFound", ...
        "LTspice was not found. Install LTspice or pass Executable=... explicitly.");
end
executable = match(1);
end
