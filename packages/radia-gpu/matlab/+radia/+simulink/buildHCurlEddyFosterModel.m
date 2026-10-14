function modelPath = buildHCurlEddyFosterModel(modelName, model, options)
%BUILDHCURLEDDYFOSTERMODEL Build a Simulink block for a Foster HCurl model.
%   The block input is u=-d(coil_current)/dt, and the output is the reduced
%   port response Bm.'*z of the diagonal Foster state z.  Use
%   makeHCurlEddyFosterModel or loadHCurlEddyFosterModel to create MODEL.

arguments
    modelName (1,1) string
    model (1,1) struct
    options.Block (1,1) string = "standard"
    options.StopTime_s (1,1) double {mustBePositive} = 1.0
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the Foster HCurl model.");
end
radia.internal.validateHCurlEddyFosterModel(model);
if ~ismember(options.Block, ["standard", "radia-mex"])
    error("radia:simulink:HCurlFosterBlock", ...
        "Block must be 'standard' or 'radia-mex'.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "close the existing model before rebuilding it.");
end

new_system(char(modelName));
root = string(modelName);
set_param(char(root), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(model.sample_time_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));

add_block("simulink/Ports & Subsystems/In1", root + "/minus_dI_dt", ...
    "Port", "1", "Position", [35 105 65 125]);
if options.Block == "standard"
    add_block("simulink/Discrete/Discrete State-Space", root + "/HCurlFoster", ...
        "A", mat2str(model.Ad, 17), "B", mat2str(model.Bd, 17), ...
        "C", mat2str(model.Cd, 17), "D", mat2str(model.Dd, 17), ...
        "X0", mat2str(model.x0, 17), ...
        "SampleTime", num2str(model.sample_time_s, 17), ...
        "Position", [145 85 380 145]);
else
    modelWorkspace = get_param(root, "ModelWorkspace");
    nativeModel = model;
    nativeModel.A = model.Ad;
    nativeModel.B = model.Bd;
    nativeModel.C = model.Cd;
    nativeModel.D = model.Dd;
    assignin(modelWorkspace, "radia_hcurl_eddy_foster_model", nativeModel);
    add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
        root + "/HCurlFoster", ...
        "FunctionName", "radia_state_space_mex_sfunction", ...
        "Parameters", "radia_hcurl_eddy_foster_model", ...
        "Position", [145 85 380 145]);
end
add_block("simulink/Ports & Subsystems/Out1", root + "/port_response", ...
    "Port", "1", "Position", [470 105 550 125]);
add_line(root, "minus_dI_dt/1", "HCurlFoster/1");
add_line(root, "HCurlFoster/1", "port_response/1");
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
