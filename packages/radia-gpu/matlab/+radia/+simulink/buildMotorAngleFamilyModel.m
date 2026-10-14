function modelPath = buildMotorAngleFamilyModel(modelName, family, options)
%BUILDMOTORANGLEFAMILYMODEL Build a native-MEX periodic motor Simulink model.

arguments
    modelName (1,1) string
    family (1,1) struct
    options.StopTime_s (1,1) double {mustBeFinite, mustBePositive} = 1.0e-2
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the periodic motor model.");
end
if ~isfield(family, "schema") || ...
        string(family.schema) ~= "radia.motor.periodic-angle-family.v1"
    error("radia:simulink:MotorAngleFamily", ...
        "family must come from makeMotorAngleFamily.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "Close the existing model before rebuilding it.");
end

new_system(char(modelName));
root = string(modelName);
set_param(char(root), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(family.sample_time_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));
workspace = get_param(char(root), "ModelWorkspace");
assignin(workspace, "radia_motor_angle_family", family);

add_block("simulink/Ports & Subsystems/In1", root + "/mechanical_angle_rad", ...
    "Port", "1", "Position", [30 65 60 85]);
add_block("simulink/Ports & Subsystems/In1", root + "/model_inputs", ...
    "Port", "2", "Position", [30 125 60 145]);
if family.input_count > 1
    set_param(root + "/model_inputs", "PortDimensions", num2str(family.input_count));
end
add_block("simulink/Signal Routing/Mux", root + "/InputMux", ...
    "Inputs", "2", "Position", [105 70 130 140]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    root + "/MotorAngleFamily", ...
    "FunctionName", "radia_motor_angle_family_mex_sfunction", ...
    "Parameters", "radia_motor_angle_family", ...
    "Position", [190 60 430 150]);
add_block("simulink/Ports & Subsystems/Out1", root + "/outputs_and_torque", ...
    "Port", "1", "Position", [490 95 600 115]);
add_line(root, "mechanical_angle_rad/1", "InputMux/1");
add_line(root, "model_inputs/1", "InputMux/2");
add_line(root, "InputMux/1", "MotorAngleFamily/1");
add_line(root, "MotorAngleFamily/1", "outputs_and_torque/1");
set_param(root, "ModelBrowserVisibility", "off");
assignin(workspace, "radia_motor_angle_family_contract", struct( ...
    "schema", family.schema, "backend", "native-mex-periodic-family", ...
    "python_per_step", false, "matlab_matrix_algebra_per_step", false));

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
