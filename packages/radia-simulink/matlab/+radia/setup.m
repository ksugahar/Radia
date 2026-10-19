function info = setup(options)
%SETUP Configure explicit Python/batch composition; never load CPU radia_mex.
arguments
    options.PythonExecutable (1,1) string = defaultPythonExecutable()
    options.RequireMex (1,1) logical = true
    options.ConfigureSimulinkFileGeneration (1,1) logical = true
    options.SimulinkFileGenerationRoot (1,1) string = ""
    options.Verbose (1,1) logical = false
    options.Force (1,1) logical = false
end
if options.RequireMex
    error("radia:mex:Retired", ...
        "CPU radia_mex has been removed. Use explicit Python CPU workflows or supported radia-simulink mexcuda operations.");
end
matlabDir=fileparts(fileparts(mfilename('fullpath')));
if ~contains(path,matlabDir), addpath(matlabDir); end
radia.internal.pythonProcessPath("capture");
radia.internal.selectedPythonExecutable("set",options.PythonExecutable);
fileGenerationInfo=struct("available",false,"changed",false,"reason","disabled");
if options.ConfigureSimulinkFileGeneration
    fileGenerationInfo=radia.simulink.configureFileGeneration( ...
        RootDirectory=options.SimulinkFileGenerationRoot,Verbose=options.Verbose);
end
% Each explicit Python operation validates its own interpreter and dependencies.
% Do not inject NGSolve/MKL DLL paths into MATLAB for a retired native gateway.
info=struct("matlab_dir",matlabDir,"mex_path","","mex_available",false, ...
    "mex_source_commit","","python_executable",options.PythonExecutable, ...
    "runtime_dirs",strings(0,1),"simulink_file_generation",fileGenerationInfo, ...
    "configured",true);
end
function executable=defaultPythonExecutable()
executable=string(getenv('RADIA_PYTHON_EXECUTABLE'));
if strlength(executable)==0, executable="python"; end
end
