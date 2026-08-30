function info = setup(options)
%SETUP Configure Radia's MATLAB runtime and locate the native MEX gateway.
%   INFO = radia.setup() is safe to call repeatedly. The first call adds the
%   package directory and NGSolve DLL directories to MATLAB's search path.

arguments
    options.PythonExecutable (1,1) string = "python"
    options.RequireMex (1,1) logical = true
    options.ConfigureSimulinkFileGeneration (1,1) logical = true
    options.SimulinkFileGenerationRoot (1,1) string = ""
    options.Verbose (1,1) logical = false
    options.Force (1,1) logical = false
end

persistent cachedInfo cachedPython
radia.internal.pythonProcessPath("capture");
packageDir = fileparts(mfilename("fullpath"));
matlabDir = fileparts(packageDir);
if ~contains(path, matlabDir)
    addpath(matlabDir);
end

fileGenerationInfo = struct("available", false, "changed", false, ...
    "root", "", "cache_folder", "", "codegen_folder", "", ...
    "reason", "disabled");
if options.ConfigureSimulinkFileGeneration
    fileGenerationInfo = radia.simulink.configureFileGeneration( ...
        RootDirectory=options.SimulinkFileGenerationRoot, ...
        Verbose=options.Verbose);
end

if ~options.Force && ~isempty(cachedInfo) && ...
        cachedPython == options.PythonExecutable && ...
        (~options.RequireMex || cachedInfo.mex_available) && ...
        (~options.RequireMex || exist("radia_mex", "file") == 3)
    cachedInfo.simulink_file_generation = fileGenerationInfo;
    info = cachedInfo;
    return
end

runtimeDirs = strings(0, 1);
excludedOpenMPRuntimeDirs = strings(0, 1);
if ispc
    python = char(options.PythonExecutable);
    code = "import os,sys,sysconfig,netgen;" + ...
        "p=sysconfig.get_path('purelib');" + ...
        "ds=(os.path.dirname(netgen.__file__)," + ...
        "os.path.join(p,'ngsolve_openblas')," + ...
        "os.path.join(sys.prefix,'bin')," + ...
        "os.path.join(sys.prefix,'Library','bin'),sys.prefix);" + ...
        "[print('RADIA_RUNTIME:'+os.path.abspath(d)) for d in ds]";
    command = sprintf('"%s" -c "%s"', python, code);
    [status, output] = system(command);
    if status ~= 0
        error("radia:setup:Python", ...
            "Could not locate the NGSolve runtime with %s:\n%s", python, output);
    end
    lines = splitlines(string(output));
    runtimeDirs = extractAfter(lines(startsWith(lines, "RADIA_RUNTIME:")), ...
        "RADIA_RUNTIME:");
    runtimeDirs = runtimeDirs(isfolder(runtimeDirs));
    foreignOpenMP = isfile(fullfile(runtimeDirs, "libiomp5md.dll"));
    excludedOpenMPRuntimeDirs = runtimeDirs(foreignOpenMP);
    runtimeDirs(foreignOpenMP) = [];
    if isempty(runtimeDirs)
        error("radia:setup:Runtime", ...
            "Python found NGSolve, but no usable runtime directories were returned.");
    end

    current = split(string(getenv("PATH")), pathsep);
    current(current == "") = [];
    managedRuntimeDirs = [runtimeDirs; excludedOpenMPRuntimeDirs];
    keep = true(size(current));
    for i = 1:numel(current)
        keep(i) = ~any(strcmpi(current(i), managedRuntimeDirs));
    end
    setenv("PATH", char(strjoin([runtimeDirs; current(keep)], pathsep)));
end

 mexPath = which("radia_mex");
 mexAvailable = exist("radia_mex", "file") == 3;
if options.RequireMex && ~mexAvailable
    error("radia:setup:MissingMex", ...
        "radia_mex is not built. Run Build.ps1 -MatlabMexOnly or configure CMake with RADIA_BUILD_MATLAB_MEX=ON.");
end

info = struct( ...
    "matlab_dir", matlabDir, ...
    "mex_path", string(mexPath), ...
    "mex_available", mexAvailable, ...
    "python_executable", options.PythonExecutable, ...
    "runtime_dirs", runtimeDirs, ...
    "excluded_openmp_runtime_dirs", excludedOpenMPRuntimeDirs, ...
    "simulink_file_generation", fileGenerationInfo, ...
    "configured", true);
if options.Verbose
    fprintf("Radia MATLAB runtime ready\n  MEX: %s\n  Python: %s\n", ...
        string(mexPath), options.PythonExecutable);
end
cachedInfo = info;
cachedPython = options.PythonExecutable;
end
