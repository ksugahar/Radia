function library = install_radia_simulink(options)
%INSTALL_RADIA_SIMULINK Configure Radia and optionally open its block library.
arguments
    options.OpenLibrary (1,1) logical = true
    options.PythonExecutable (1,1) string = "python"
end
matlabRoot = fileparts(mfilename("fullpath"));
addpath(matlabRoot, "-begin");
radia.setup(PythonExecutable=options.PythonExecutable, ...
    RequireMex=true,Force=true);
library = fullfile(matlabRoot, "radia_simulink_library.slx");
if options.OpenLibrary
    load_system(library);
    open_system("radia_simulink_library");
end
end
