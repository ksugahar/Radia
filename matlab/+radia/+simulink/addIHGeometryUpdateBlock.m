function blockPath = addIHGeometryUpdateBlock(modelName, options)
%ADDIHGEOMETRYUPDATEBLOCK Place the IH geometry watch/rebuild block.
%   blockPath = addIHGeometryUpdateBlock(modelName) adds a masked
%   "Geometry Update" block and installs a model InitFcn hook that runs
%   radia.simulink.updateIHGeometry on every diagram update /
%   simulation start.  The mask holds the two geometry inputs
%   (workpiece .vol and coil .step/.vol -- roles are classified by
%   extension, so a crossed pair is repaired), the assemble command
%   that rebuilds the operators, the configuration file that command
%   writes, and an auto-rebuild switch; a "Rebuild now" button forces a
%   rebuild regardless of the fingerprint. Browse buttons beside both
%   geometry fields select and assign absolute paths without rebuilding
%   until the next diagram update, simulation start, or explicit button.
%
%   With the block configured, re-pointing or editing a geometry file
%   is enough: the next update detects the content change through the
%   SHA-256 fingerprint sidecar, reruns the assemble command at the
%   explicit-update boundary (never per time step), and reloads the
%   refreshed configuration through configureIHNativeModel.

arguments
    modelName (1,1) string
    options.BlockName (1,1) string = "Geometry Update"
    options.Position (1,4) double = [185 345 365 425]
end

blockPath = modelName + "/" + options.BlockName;
add_block("simulink/Ports & Subsystems/Subsystem", blockPath, ...
    Position=options.Position);
delete_line(blockPath, "In1/1", "Out1/1");
delete_block(blockPath + "/In1");
delete_block(blockPath + "/Out1");
set_param(blockPath, "Tag", "RadiaIHGeometryUpdate");

mask = Simulink.Mask.create(blockPath);
mask.Description = "Watches the workpiece/coil geometry files and " + ...
    "rebuilds the IH operators when their content changes.  Set " + ...
    "EITHER the assemble function (a MATLAB .m on the path, called " + ...
    "as fcn(wpVol, coilFile, configFile) -- preferred, in-process, " + ...
    "errors propagate) OR the assemble command (shell, for CLI " + ...
    "assemblers); it must write config_file (MAT/JSON, native IH " + ...
    "configuration) and runs only at diagram update / simulation " + ...
    "start or through the Rebuild now button, never per time step.  " + ...
    "Use absolute file paths.  Leave every field empty to keep the " + ...
    "block inert.";
wpParameter = mask.addParameter(Type="edit", Name="wp_vol", ...
    Prompt="Workpiece mesh (.vol / .vol.gz)", Value="", Evaluate="off");
wpParameter.ShowTooltip = "on";
wpBrowse = mask.addDialogControl("pushbutton", "browse_wp_vol");
wpBrowse.Prompt = "Browse...";
wpBrowse.Tooltip = "Select the workpiece Netgen .vol or .vol.gz mesh.";
wpBrowse.Callback = ...
    "radia.simulink.browseIHGeometryFile(string(gcb), 'workpiece');";
wpBrowse.Row = "current";
wpBrowse.HorizontalStretch = "off";

coilParameter = mask.addParameter(Type="edit", Name="coil_file", ...
    Prompt="Coil geometry (.step/.stp PEEC / .vol BEM-A)", Value="", ...
    Evaluate="off");
coilParameter.ShowTooltip = "on";
coilBrowse = mask.addDialogControl("pushbutton", "browse_coil_file");
coilBrowse.Prompt = "Browse...";
coilBrowse.Tooltip = ...
    "Select STEP CAD for PEEC or a .vol/.vol.gz mesh for BEM-A.";
coilBrowse.Callback = ...
    "radia.simulink.browseIHGeometryFile(string(gcb), 'coil');";
coilBrowse.Row = "current";
coilBrowse.HorizontalStretch = "off";

mask.addParameter(Type="edit", Name="assemble_fcn", ...
    Prompt="Assemble function (MATLAB, preferred): " + ...
    "fcn(wpVol, coilFile, configFile)", Value="", Evaluate="off");
mask.addParameter(Type="edit", Name="assemble_command", ...
    Prompt="Assemble command (shell, alternative; writes config_file)", ...
    Value="", Evaluate="off");
mask.addParameter(Type="edit", Name="config_file", ...
    Prompt="Configuration MAT/JSON the command writes", Value="", ...
    Evaluate="off");
mask.addParameter(Type="checkbox", Name="auto_rebuild", ...
    Prompt="Rebuild automatically when files change", Value="on");
button = mask.addDialogControl("pushbutton", "rebuild_now");
button.Prompt = "Rebuild now";
button.Callback = ...
    "radia.simulink.updateIHGeometry(string(bdroot(gcb)), Force=true);";
mask.Display = "disp('Geometry Update');";

set_param(blockPath, "UserData", struct( ...
    "role", "geometry-watch-rebuild", ...
    "fingerprint", "sha256-content+command sidecar", ...
    "execution", "model-init-or-explicit-button", ...
    "python_per_step", false), ...
    "UserDataPersistent", "on");

hook = "radia.simulink.updateIHGeometry(string(bdroot));";
current = string(get_param(modelName, "InitFcn"));
if ~contains(current, "radia.simulink.updateIHGeometry")
    if strlength(strtrim(current)) == 0
        set_param(modelName, "InitFcn", hook);
    else
        set_param(modelName, "InitFcn", current + newline + hook);
    end
end
end
