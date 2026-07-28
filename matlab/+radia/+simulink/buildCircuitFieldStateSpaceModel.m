function modelPath = buildCircuitFieldStateSpaceModel(modelName, model, options)
%BUILDCIRCUITFIELDSTATESPACEMODEL Build a native-MEX circuit Simulink model.

arguments
    modelName (1,1) string
    model (1,1) struct
    options.StopTime_s (1,1) double {mustBePositive} = 1.0e-2
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the circuit-field model.");
end
if ~isfield(model, "schema") || ...
        string(model.schema) ~= "radia.circuit-field.state-space.v1"
    error("radia:simulink:CircuitFieldModel", ...
        "model must come from makeCircuitFieldStateSpace.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "Close the existing model before rebuilding it.");
end

new_system(char(modelName));
root = string(modelName);
set_param(char(root), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(model.sample_time_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));

add_block("simulink/Ports & Subsystems/In1", root + "/voltage", ...
    "Port", "1", "Position", [35 105 65 125]);
modelWorkspace = get_param(root, "ModelWorkspace");
nativeModel = model;
nativeModel.A = model.Ad;
nativeModel.B = model.Bd;
nativeModel.C = model.Cd;
nativeModel.D = model.Dd;
assignin(modelWorkspace, "radia_circuit_field_model", nativeModel);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    root + "/Circuit Field ROM", ...
    "FunctionName", "radia_state_space_mex_sfunction", ...
    "Parameters", "radia_circuit_field_model", ...
    "Position", [145 75 385 155]);
add_block("simulink/Ports & Subsystems/Out1", root + "/current_and_flux", ...
    "Port", "1", "Position", [475 105 575 125]);
add_line(root, "voltage/1", "Circuit Field ROM/1");
add_line(root, "Circuit Field ROM/1", "current_and_flux/1");
set_param(root, "ModelBrowserVisibility", "off");
assignin(modelWorkspace, "radia_circuit_field_contract", struct( ...
    "schema", model.schema, "backend", "native-mex-sfunction", ...
    "python_per_step", false, "field_factorization_per_step", false));

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
