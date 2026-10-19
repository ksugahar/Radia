function executable = resolveAssemblyPython(explicit)
%RESOLVEASSEMBLYPYTHON Select Python for an explicit offline assembly.
%   A nonblank caller override is returned unchanged. Otherwise reuse the
%   interpreter validated by radia.setup, configuring the default runtime only
%   when setup has not yet run in this MATLAB process.

arguments
    explicit (1,1) string = ""
end
if strlength(strtrim(explicit)) > 0
    executable = explicit;
    return
end
executable = radia.internal.selectedPythonExecutable("get");
if strlength(executable) == 0
    setupInfo = radia.setup(RequireMex=false, ...
        ConfigureSimulinkFileGeneration=false);
    executable = setupInfo.python_executable;
end
end
