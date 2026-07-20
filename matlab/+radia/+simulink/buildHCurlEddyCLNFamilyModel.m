function modelPath = buildHCurlEddyCLNFamilyModel(modelName, family, options)
%BUILDHCURLEDDYCLNFAMILYMODEL Build a moving common-basis HCurl CLN model.
%   Inputs are minus_dI_dt, height_offset_m, and coil_current_A.  The single
%   output vector contains port response followed by the three force
%   components.  The family must be created by loadHCurlEddyCLNFamily or
%   contain the same validated schema.

arguments
    modelName (1,1) string
    family (1,1) struct
    options.StopTime_s (1,1) double {mustBePositive} = 1.0
    options.Save (1,1) logical = true
    options.Open (1,1) logical = false
end

if exist("new_system", "file") ~= 2 && exist("new_system", "builtin") ~= 5
    error("radia:simulink:MissingSimulink", ...
        "Simulink is required to build the HCurl CLN family model.");
end
if ~isfield(family, "schema") || family.schema ~= ...
        "radia.hcurl.eddy_cln.family.v1"
    error("radia:simulink:HCurlCLNFamily", ...
        "family must come from loadHCurlEddyCLNFamily.");
end
if bdIsLoaded(char(modelName))
    error("radia:simulink:ModelLoaded", ...
        "close the existing model before rebuilding it.");
end

new_system(char(modelName));
root = string(modelName);
set_param(char(root), "Solver", "FixedStepDiscrete", ...
    "FixedStep", num2str(family.sample_time_s, 17), ...
    "StopTime", num2str(options.StopTime_s, 17));
workspace = get_param(char(root), "ModelWorkspace");
assignin(workspace, "radia_hcurl_eddy_cln_family", family);

nPort = family.port_count;
nInput = 2 * nPort + 1;
add_block("simulink/Ports & Subsystems/In1", root + "/minus_dI_dt", ...
    "Port", "1", "Position", [30 70 60 90]);
add_block("simulink/Ports & Subsystems/In1", root + "/height_offset_m", ...
    "Port", "2", "Position", [30 120 60 140]);
add_block("simulink/Ports & Subsystems/In1", root + "/coil_current_A", ...
    "Port", "3", "Position", [30 170 60 190]);
if nPort > 1
    set_param(root + "/minus_dI_dt", "PortDimensions", num2str(nPort));
    set_param(root + "/coil_current_A", "PortDimensions", num2str(nPort));
end
add_block("simulink/Signal Routing/Mux", root + "/InputMux", ...
    "Inputs", "3", "Position", [115 85 140 175]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    root + "/HCurlCLNFamily", ...
    "FunctionName", "radia_hcurl_eddy_cln_family_sfunction", ...
    "Parameters", "radia_hcurl_eddy_cln_family", ...
    "Position", [200 75 440 185]);
add_block("simulink/Ports & Subsystems/Out1", root + "/response_and_force", ...
    "Port", "1", "Position", [500 115 620 135]);
add_line(root, "minus_dI_dt/1", "InputMux/1");
add_line(root, "height_offset_m/1", "InputMux/2");
add_line(root, "coil_current_A/1", "InputMux/3");
add_line(root, "InputMux/1", "HCurlCLNFamily/1");
add_line(root, "HCurlCLNFamily/1", "response_and_force/1");
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
