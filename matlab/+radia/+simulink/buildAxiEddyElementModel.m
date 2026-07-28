function modelPath = buildAxiEddyElementModel(modelName, model, options)
%BUILDAXIEDDYELEMENTMODEL Build a native-MEX Simulink Q2 eddy model.

arguments
    modelName (1,1) string
    model (1,1) struct
    options.StopTime_s (1,1) double {mustBePositive} = 1.0e-3
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the axifem eddy model.");
end
if ~isfield(model, "schema") || ...
        string(model.schema) ~= "radia.axifem.q2_eddy.state_space.v1"
    error("radia:simulink:AxiEddyModel", ...
        "model must come from makeAxiEddyElementModel.");
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

add_block("simulink/Ports & Subsystems/In1", root + "/current", ...
    "Port", "1", "Position", [35 105 65 125]);
modelWorkspace = get_param(root, "ModelWorkspace");
nativeModel = model;
nativeModel.A = model.Ad;
nativeModel.B = model.Bd;
nativeModel.C = model.Cd;
nativeModel.D = model.Dd;
assignin(modelWorkspace, "radia_axifem_eddy_model", nativeModel);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    root + "/Q2 Eddy Element", ...
    "FunctionName", "radia_state_space_mex_sfunction", ...
    "Parameters", "radia_axifem_eddy_model", ...
    "Position", [145 75 385 155]);
add_block("simulink/Ports & Subsystems/Out1", root + "/Aphi_nodes", ...
    "Port", "1", "Position", [475 105 555 125]);
add_line(root, "current/1", "Q2 Eddy Element/1");
add_line(root, "Q2 Eddy Element/1", "Aphi_nodes/1");
set_param(root, "ModelBrowserVisibility", "off");
assignin(modelWorkspace, "radia_axifem_eddy_contract", struct( ...
    "schema", model.schema, "backend", "native-mex-sfunction", ...
    "python_per_step", false, "element_order", "Q2"));

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
