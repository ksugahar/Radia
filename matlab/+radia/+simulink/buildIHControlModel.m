function modelPath = buildIHControlModel(modelName, plant, options)
%BUILD IHCONTROLMODEL Create a small Simulink model around an IH plant.
%   radia.simulink.buildIHControlModel(NAME,PLANT) creates a fixed-step
%   model with inputs power_W and ambient_temperature_K and outputs matching
%   PLANT.output_names. Set IncludePID=true to add a temperature-setpoint
%   input and a bounded discrete PID controller in front of the plant. The
%   electromagnetic power provider remains external: it may be a Radia LUT,
%   VIM/FEM co-simulation block, or a measured drive signal. Set
%   PlantBlock="radia-sfunction" to use the Radia-owned Level-2 MATLAB
%   S-function boundary instead of the standard Discrete State-Space block.

arguments
    modelName (1,1) string
    plant (1,1) struct
    options.IncludePID (1,1) logical = false
    options.Kp (1,1) double {mustBeFinite} = 1.0
    options.Ki (1,1) double {mustBeFinite} = 0.0
    options.Kd (1,1) double {mustBeFinite} = 0.0
    options.PowerLowerBound_W (1,1) double {mustBeNonnegative} = 0.0
    options.PowerUpperBound_W (1,1) double {mustBePositive} = 1.0e6
    options.StopTime_s (1,1) double {mustBePositive} = 10.0
    options.PlantBlock (1,1) string = "standard"
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the IH control model.");
end
if ~isfield(plant, "schema") || plant.schema ~= "radia.ih.simulink.plant.v1"
    error("radia:simulink:InvalidPlant", "plant is not a Radia IH plant.");
end
if ~ismember(options.PlantBlock, ["standard", "radia-sfunction", "radia-mex"])
    error("radia:simulink:PlantBlock", ...
        "PlantBlock must be 'standard', 'radia-sfunction', or 'radia-mex'.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "close the existing model before rebuilding it.");
end

new_system(char(modelName));
set_param(char(modelName), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(plant.sample_time_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));

root = string(modelName);
add_block("simulink/Signal Routing/Mux", root + "/InputMux", ...
    "Inputs", "2", "Position", [535 95 560 165]);
if options.PlantBlock == "standard"
    add_block("simulink/Discrete/Discrete State-Space", root + "/IHPlant", ...
        "A", mat2str(plant.A, 17), "B", mat2str(plant.B, 17), ...
        "C", mat2str(plant.C, 17), "D", mat2str(plant.D, 17), ...
        "X0", mat2str(plant.x0, 17), ...
        "SampleTime", num2str(plant.sample_time_s, 17), ...
        "Position", [620 80 800 170]);
elseif options.PlantBlock == "radia-sfunction"
    modelWorkspace = get_param(root, "ModelWorkspace");
    assignin(modelWorkspace, "radia_ih_plant", plant);
    add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
        root + "/IHPlant", ...
        "FunctionName", "radia_ih_plant_sfunction", ...
        "Parameters", "radia_ih_plant", ...
        "Position", [620 80 800 170]);
else
    modelWorkspace = get_param(root, "ModelWorkspace");
    assignin(modelWorkspace, "radia_ih_plant", plant);
    add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
        root + "/IHPlant", ...
        "FunctionName", "radia_state_space_mex_sfunction", ...
        "Parameters", "radia_ih_plant", ...
        "Position", [620 80 800 170]);
end
add_block("simulink/Signal Routing/Demux", root + "/OutputDemux", ...
    "Outputs", num2str(size(plant.C, 1)), "Position", [855 60 880 190]);

if options.IncludePID
    add_block("simulink/Math Operations/Sum", root + "/TemperatureError", ...
        "Inputs", "+-", "Position", [130 70 160 110]);
    add_block("simulink/Discrete/Discrete PID Controller", root + "/TemperaturePID", ...
        "P", num2str(options.Kp, 17), "I", num2str(options.Ki, 17), ...
        "D", num2str(options.Kd, 17), ...
        "SampleTime", num2str(plant.sample_time_s, 17), ...
        "Position", [210 60 340 120]);
    add_block("simulink/Discontinuities/Saturation", root + "/PowerSaturation", ...
        "UpperLimit", num2str(options.PowerUpperBound_W, 17), ...
        "LowerLimit", num2str(options.PowerLowerBound_W, 17), ...
        "Position", [390 70 480 110]);
    add_block("simulink/Discrete/Unit Delay", root + "/TemperatureFeedbackDelay", ...
        "InitialCondition", num2str(plant.x0(1), 17), ...
        "SampleTime", num2str(plant.sample_time_s, 17), ...
        "Position", [350 260 420 300]);
    add_block("simulink/Ports & Subsystems/In1", root + "/temperature_setpoint_K", ...
        "Port", "1", "Position", [25 80 55 100]);
    add_block("simulink/Ports & Subsystems/In1", root + "/ambient_temperature_K", ...
        "Port", "2", "Position", [425 165 455 185]);
    add_line(root, "temperature_setpoint_K/1", "TemperatureError/1");
    add_line(root, "TemperatureError/1", "TemperaturePID/1");
    add_line(root, "TemperaturePID/1", "PowerSaturation/1");
    add_line(root, "PowerSaturation/1", "InputMux/1");
else
    add_block("simulink/Ports & Subsystems/In1", root + "/power_W", ...
        "Port", "1", "Position", [425 95 455 115]);
    add_block("simulink/Ports & Subsystems/In1", root + "/ambient_temperature_K", ...
        "Port", "2", "Position", [425 145 455 165]);
    add_line(root, "power_W/1", "InputMux/1");
end

add_line(root, "ambient_temperature_K/1", "InputMux/2");
add_line(root, "InputMux/1", "IHPlant/1");
add_line(root, "IHPlant/1", "OutputDemux/1");
for k = 1:size(plant.C, 1)
    name = char(plant.output_names(k));
    blockName = matlab.lang.makeValidName(name);
    add_block("simulink/Ports & Subsystems/Out1", root + "/" + blockName, ...
        "Port", num2str(k), "Position", [980 35 + 50 * k 1060 55 + 50 * k]);
    add_line(root, "OutputDemux/" + k, blockName + "/1");
    set_param(root + "/" + blockName, "Name", name);
end

if options.IncludePID
    add_line(root, "OutputDemux/1", "TemperatureFeedbackDelay/1", ...
        "autorouting", "smart");
    add_line(root, "TemperatureFeedbackDelay/1", "TemperatureError/2", ...
        "autorouting", "smart");
end

set_param(root, "ModelBrowserVisibility", "off");
if options.Open
    open_system(root);
end
if options.Save
    save_system(root);
    modelPath = string(get_param(root, "FileName"));
else
    modelPath = "";
end
end
