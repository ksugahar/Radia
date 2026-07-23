function modelPath = buildIHNativeModel(options)
%BUILDIHNATIVEMODEL Build the standalone native Eddy/Thermal IH preview.
arguments
    options.ModelName (1,1) string = "radia_ih"
    options.OutputDirectory (1,1) string = ""
    options.Config (1,1) struct = struct()
    options.Open (1,1) logical = false
end
radia.simulink.requireIHNativeRuntime();
matlabRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
if strlength(options.OutputDirectory) == 0
    options.OutputDirectory = matlabRoot;
end
if ~isfolder(options.OutputDirectory)
    mkdir(options.OutputDirectory);
end
if bdIsLoaded(options.ModelName)
    close_system(options.ModelName, 0);
end
config = options.Config;
if isempty(fieldnames(config))
    config = radia.simulink.makeIHNativeSmokeConfig();
end

new_system(options.ModelName);
workspace = get_param(options.ModelName, "ModelWorkspace");
workspace.assignin("radia_ih_config", config);

add_block("simulink/Sources/Constant", ...
    options.ModelName + "/Coil Current", ...
    Value="1", Position=[35 55 95 85]);
add_block("simulink/Sources/Ramp", ...
    options.ModelName + "/Workpiece Angle", ...
    slope="0", start="0", InitialOutput="0", ...
    Position=[35 135 95 165]);
add_block("simulink/Sources/Constant", ...
    options.ModelName + "/Ambient Temperature", ...
    Value="293.15", Position=[35 235 95 265]);

add_block("simulink/User-Defined Functions/S-Function", ...
    options.ModelName + "/Eddy", ...
    FunctionName="radia_ih_eddy_sfun", ...
    Parameters="radia_ih_config", Position=[195 55 365 125]);
add_block("simulink/User-Defined Functions/S-Function", ...
    options.ModelName + "/Thermal", ...
    FunctionName="radia_ih_thermal_sfun", ...
    Parameters="radia_ih_config", Position=[470 145 640 215]);

add_block("simulink/Ports & Subsystems/Out1", ...
    options.ModelName + "/heat_density_W_per_m3", ...
    Port="1", Position=[720 70 750 90]);
add_block("simulink/Ports & Subsystems/Out1", ...
    options.ModelName + "/temperature_K", ...
    Port="2", Position=[720 165 750 185]);
add_block("simulink/Sinks/Scope", ...
    options.ModelName + "/Heat Density", Position=[805 55 900 105]);
add_block("simulink/Sinks/Scope", ...
    options.ModelName + "/Temperature", Position=[805 145 900 195]);

parameterPath = options.ModelName + "/IH Parameters";
add_block("simulink/Ports & Subsystems/Subsystem", parameterPath, ...
    Position=[185 255 365 315]);
delete_line(parameterPath, "In1/1", "Out1/1");
delete_block(parameterPath + "/In1");
delete_block(parameterPath + "/Out1");
mask = Simulink.Mask.create(parameterPath);
mask.Description = "Native IH preview configuration with preassembled operators. MAT files contain config or radia_ih_config; JSON is also supported.";
configParameter = mask.addParameter(Type="edit", Name="config_file", ...
    Prompt="IH configuration MAT/JSON", Value="''", Evaluate="on");
configParameter.Callback = ...
    "radia.simulink.configureIHNativeModel(string(bdroot(gcb)),string(config_file));";
mask.Display = "disp('IH Parameters');";
set_param(parameterPath, "UserData", struct( ...
    "backend", "native-mex-sfunction", ...
    "release_channel", "preview", ...
    "operator_assembly", "preassembled", ...
    "python_fallback", false, ...
    "configuration", "model-workspace-or-file"), ...
    "UserDataPersistent", "on");

connect(options.ModelName, "Coil Current/1", "Eddy/1");
connect(options.ModelName, "Workpiece Angle/1", "Eddy/2");
connect(options.ModelName, "Workpiece Angle/1", "Thermal/3");
connect(options.ModelName, "Ambient Temperature/1", "Thermal/2");
connect(options.ModelName, "Eddy/1", "Thermal/1");
connect(options.ModelName, "Eddy/1", "heat_density_W_per_m3/1");
connect(options.ModelName, "Eddy/1", "Heat Density/1");
connect(options.ModelName, "Thermal/1", "Eddy/3");
connect(options.ModelName, "Thermal/1", "temperature_K/1");
connect(options.ModelName, "Thermal/1", "Temperature/1");

set_param(options.ModelName, "StopTime", "1", ...
    "SaveOutput", "on", "OutputSaveName", "yout");
radia.simulink.configureIHNativeModel(options.ModelName);
modelPath = fullfile(options.OutputDirectory, options.ModelName + ".slx");
save_system(options.ModelName, modelPath);
if options.Open
    open_system(options.ModelName);
else
    close_system(options.ModelName, 0);
end

function connect(modelName, source, destination)
add_line(modelName, source, destination, "autorouting", "smart");
end
end
