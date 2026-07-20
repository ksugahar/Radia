function modelPath = buildTeam28CLNModel(modelName, lut, options)
%BUILDTEAM28CLNMODEL Build a Simulink model around the TEAM 28 CLN LUT.
%   Inputs are height_offset_m and coil_current_A.  Outputs are physical
%   signed force_N, positive upward_lift_N, and force slope.  The model is a
%   control-oriented 50 Hz benchmark model; it does not extrapolate to a new
%   frequency or replace the high-fidelity HCurl-VIM gate.

arguments
    modelName (1,1) string
    lut (1,1) struct
    options.SampleTime_s (1,1) double {mustBePositive} = 1.0e-3
    options.StopTime_s (1,1) double {mustBePositive} = 1.0
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the TEAM 28 model.");
end
if ~isfield(lut, "schema") || lut.schema ~= "radia.team28.cln_lut.v1"
    error("radia:simulink:Team28LUT", "lut must come from makeTeam28CLNLUT.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "close the existing model before rebuilding it.");
end

new_system(char(modelName));
root = string(modelName);
set_param(char(root), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(options.SampleTime_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));

modelWorkspace = get_param(char(root), "ModelWorkspace");
lut.sample_time_s = options.SampleTime_s;
assignin(modelWorkspace, "radia_team28_cln_lut", lut);

add_block("simulink/Ports & Subsystems/In1", root + "/height_offset_m", ...
    "Port", "1", "Position", [30 75 60 95]);
add_block("simulink/Ports & Subsystems/In1", root + "/coil_current_A", ...
    "Port", "2", "Position", [30 135 60 155]);
add_block("simulink/Signal Routing/Mux", root + "/InputMux", ...
    "Inputs", "2", "Position", [135 90 160 150]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    root + "/Team28CLNLUT", ...
    "FunctionName", "radia_team28_cln_lut_sfunction", ...
    "Parameters", "radia_team28_cln_lut", ...
    "Position", [220 75 430 165]);
add_block("simulink/Signal Routing/Demux", root + "/OutputDemux", ...
    "Outputs", "3", "Position", [485 60 510 180]);

outputs = ["force_N", "upward_lift_N", "force_slope_N_per_m"];
add_line(root, "height_offset_m/1", "InputMux/1");
add_line(root, "coil_current_A/1", "InputMux/2");
add_line(root, "InputMux/1", "Team28CLNLUT/1");
add_line(root, "Team28CLNLUT/1", "OutputDemux/1");
for k = 1:numel(outputs)
    name = outputs(k);
    blockName = matlab.lang.makeValidName(name);
    add_block("simulink/Ports & Subsystems/Out1", root + "/" + blockName, ...
        "Port", num2str(k), "Position", [595 45 + 55 * k 700 65 + 55 * k]);
    add_line(root, "OutputDemux/" + k, blockName + "/1");
    set_param(root + "/" + blockName, "Name", name);
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
