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
workspace.assignin("radia_ih_geometry_revision", 0);
radia.simulink.makeIHMonitorBusObject();

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

add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    options.ModelName + "/Eddy", ...
    FunctionName="radia_ih_eddy_sfun", ...
    Parameters="radia_ih_config", Position=[195 55 365 125]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    options.ModelName + "/Thermal", ...
    FunctionName="radia_ih_thermal_sfun", ...
    Parameters="radia_ih_config", Position=[470 145 640 215]);

add_block("simulink/Ports & Subsystems/Out1", ...
    options.ModelName + "/heat_density_W_per_m3", ...
    Port="1", Position=[720 70 750 90]);
add_block("simulink/Ports & Subsystems/Out1", ...
    options.ModelName + "/temperature_K", ...
    Port="2", Position=[720 165 750 185]);

parameterPath = options.ModelName + "/IH Parameters";
add_block("simulink/Ports & Subsystems/Subsystem", parameterPath, ...
    Position=[185 255 365 315]);
delete_line(parameterPath, "In1/1", "Out1/1");
delete_block(parameterPath + "/In1");
delete_block(parameterPath + "/Out1");
mask = Simulink.Mask.create(parameterPath);
mask.Description = "Native IH preview configuration with preassembled operators. MAT files contain config or radia_ih_config; JSON is also supported.";
configParameter = mask.addParameter(Type="edit", Name="config_file", ...
    Prompt="IH configuration MAT/JSON", Value="", Evaluate="off");
% Thin callback: the .slx stores only a named .m entry point (inline
% mask-callback code is invisible to diffs and cannot be tested; and
% mask callbacks run in the base workspace where bare mask variables
% like config_file do not exist).
configParameter.Callback = "radia.simulink.onIHConfigFileChanged(gcb);";
mask.Display = "disp('IH Parameters');";
set_param(parameterPath, "UserData", struct( ...
    "backend", "matlab-level2+radia-mex-handles", ...
    "release_channel", "preview", ...
    "operator_assembly", "preassembled", ...
    "python_fallback", false, ...
    "configuration", "model-workspace-or-file"), ...
    "UserDataPersistent", "on");

% Geometry watch/rebuild: re-pointing or editing the .vol/.step inputs
% triggers the assemble command at the next update (InitFcn hook).
radia.simulink.addIHGeometryUpdateBlock(options.ModelName, ...
    Position=[185 345 365 425]);

connect(options.ModelName, "Coil Current/1", "Eddy/1");
connect(options.ModelName, "Workpiece Angle/1", "Eddy/2");
connect(options.ModelName, "Workpiece Angle/1", "Thermal/3");
connect(options.ModelName, "Ambient Temperature/1", "Thermal/2");
connect(options.ModelName, "Eddy/1", "Thermal/1");
connect(options.ModelName, "Eddy/1", "heat_density_W_per_m3/1");
connect(options.ModelName, "Thermal/1", "Eddy/3");
connect(options.ModelName, "Thermal/1", "temperature_K/1");
radia.simulink.addIHMonitorBus(options.ModelName, ...
    HeatSource="Eddy/1", TemperatureSource="Thermal/1", ...
    CurrentSource="Coil Current/1", AngleSource="Workpiece Angle/1", ...
    AmbientSource="Ambient Temperature/1", OutportNumber=3, ...
    Position=[700 250 850 390], OutportPosition=[930 305 960 325]);
radia.simulink.addIHMonitorDashboard(options.ModelName, ...
    Position=[1010 250 1160 330]);
connect(options.ModelName, "IH Monitor/1", "IH Dashboard/1");

set_param(options.ModelName, "StopTime", "1", ...
    "SaveOutput", "on", "OutputSaveName", "yout", ...
    "PreLoadFcn", "radia.simulink.makeIHMonitorBusObject();");
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
