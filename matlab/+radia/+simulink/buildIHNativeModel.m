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
add_block("simulink/Sinks/Scope", ...
    options.ModelName + "/Heat Density", Position=[805 55 900 105]);
add_block("simulink/Sinks/Scope", ...
    options.ModelName + "/Temperature", Position=[805 145 900 195]);
% Scope blocks default to Manual axes with YLimits [-10 10]; both IH
% signals live far outside that window (temperature ~293 K, real heat
% densities ~1e6 W/m^3), so the trace renders off-screen and the scope
% looks empty.  Follow the data instead.
scopeNames = ["Heat Density", "Temperature"];
for scopeIndex = 1:numel(scopeNames)
    scopeConfiguration = get_param( ...
        options.ModelName + "/" + scopeNames(scopeIndex), ...
        "ScopeConfiguration");
    scopeConfiguration.AxesScaling = "Auto";
end
% A real configuration carries the full field vector (measured
% 2026-08-04: 3122 temperature DOFs); a scope fed that raw vector draws
% thousands of overlapping lines and is unreadable, so each scope shows
% the [min mean max] reduction instead.  The outports still carry the
% untouched full vectors.
addVectorStatsSubsystem(options.ModelName, "Heat Stats", ...
    [700 55 770 105]);
addVectorStatsSubsystem(options.ModelName, "Temperature Stats", ...
    [700 145 770 195]);

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
% Mask callbacks run in the base workspace, where mask parameter
% variables are NOT defined -- the dialog value must be read back with
% get_param(gcb, ...).  Referencing bare config_file made every OK
% press fail with "'config_file' is not recognized".
configParameter.Callback = ...
    "radia.simulink.configureIHNativeModel(string(bdroot(gcb)), " + ...
    "string(get_param(gcb, 'config_file')));";
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
connect(options.ModelName, "Eddy/1", "Heat Stats/1");
connect(options.ModelName, "Heat Stats/1", "Heat Density/1");
connect(options.ModelName, "Thermal/1", "Eddy/3");
connect(options.ModelName, "Thermal/1", "temperature_K/1");
connect(options.ModelName, "Thermal/1", "Temperature Stats/1");
connect(options.ModelName, "Temperature Stats/1", "Temperature/1");

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

function addVectorStatsSubsystem(modelName, blockName, position)
% [min mean max] display reduction for an N-wide field vector.  The
% mean is the arithmetic mean over DOFs (a display aid; the reported
% volume-weighted T_mean stays owned by the result artifacts).
path = modelName + "/" + blockName;
add_block("simulink/Ports & Subsystems/Subsystem", path, ...
    Position=position);
delete_line(path, "In1/1", "Out1/1");
set_param(path + "/In1", "Position", [40 128 70 142]);
set_param(path + "/Out1", "Position", [420 138 450 152]);
add_block("simulink/Math Operations/MinMax", path + "/Min", ...
    Function="min", Inputs="1", Position=[170 40 210 70]);
add_block("simulink/Math Operations/Sum", path + "/Sum", ...
    Inputs="+", CollapseMode="All dimensions", ...
    Position=[170 100 210 130]);
add_block("simulink/Signal Attributes/Width", path + "/Width", ...
    Position=[170 160 210 190]);
add_block("simulink/Math Operations/Divide", path + "/Mean", ...
    Inputs="*/", Position=[250 118 290 152]);
add_block("simulink/Math Operations/MinMax", path + "/Max", ...
    Function="max", Inputs="1", Position=[170 220 210 250]);
add_block("simulink/Signal Routing/Mux", path + "/Mux", ...
    Inputs="3", Position=[340 40 350 250]);
add_line(path, "In1/1", "Min/1", "autorouting", "smart");
add_line(path, "In1/1", "Sum/1", "autorouting", "smart");
add_line(path, "In1/1", "Width/1", "autorouting", "smart");
add_line(path, "In1/1", "Max/1", "autorouting", "smart");
add_line(path, "Sum/1", "Mean/1", "autorouting", "smart");
add_line(path, "Width/1", "Mean/2", "autorouting", "smart");
set_param(add_line(path, "Min/1", "Mux/1", "autorouting", "smart"), ...
    "Name", "min");
set_param(add_line(path, "Mean/1", "Mux/2", "autorouting", "smart"), ...
    "Name", "mean");
set_param(add_line(path, "Max/1", "Mux/3", "autorouting", "smart"), ...
    "Name", "max");
add_line(path, "Mux/1", "Out1/1", "autorouting", "smart");
end
end
