function selectedPath = browseIHGeometryFile(blockPath, role, options)
%BROWSEIHGEOMETRYFILE Select and assign one IH geometry input file.
%   selectedPath = browseIHGeometryFile(blockPath, role) opens a file
%   chooser for the Geometry Update block. role is "workpiece" or
%   "coil". Workpiece accepts Netgen .vol/.vol.gz meshes; coil accepts
%   .step/.stp CAD for PEEC or .vol/.vol.gz meshes for BEM-A.
%
%   Cancelling leaves the mask parameter unchanged and returns "".
%   DialogFcn is injectable so selection behavior can be tested without
%   opening desktop UI.

arguments
    blockPath (1,1) string
    role (1,1) string
    options.DialogFcn (1,1) function_handle = @uigetfile
end

if getSimulinkBlockHandle(blockPath) < 0
    error("radia:simulink:IHGeometryBrowseBlock", ...
        "Geometry Update block does not exist: %s", blockPath);
end

role = lower(strtrim(role));
switch role
    case "workpiece"
        parameterName = "wp_vol";
        titleText = "Select workpiece mesh";
        filterSpec = { ...
            '*.vol;*.vol.gz', 'Netgen mesh (*.vol, *.vol.gz)'; ...
            '*.*', 'All files (*.*)'};
        extensions = [".vol", ".vol.gz"];
    case "coil"
        parameterName = "coil_file";
        titleText = "Select coil geometry";
        filterSpec = { ...
            '*.step;*.stp;*.vol;*.vol.gz', ...
            'Coil geometry (*.step, *.stp, *.vol, *.vol.gz)'; ...
            '*.*', 'All files (*.*)'};
        extensions = [".step", ".stp", ".vol", ".vol.gz"];
    otherwise
        error("radia:simulink:IHGeometryBrowseRole", ...
            "Geometry role must be 'workpiece' or 'coil'; got: %s", role);
end

currentValue = strtrim(string(get_param(blockPath, parameterName)));
startLocation = string(pwd);
if isfile(currentValue)
    startLocation = currentValue;
elseif strlength(currentValue) > 0
    parent = string(fileparts(currentValue));
    if isfolder(parent)
        startLocation = parent;
    end
end

[fileName, folderName] = options.DialogFcn( ...
    filterSpec, char(titleText), char(startLocation));
if isequal(fileName, 0) || isequal(folderName, 0)
    selectedPath = "";
    return
end

selectedPath = string(fullfile(string(folderName), string(fileName)));
if ~any(endsWith(lower(selectedPath), extensions))
    error("radia:simulink:IHGeometryBrowseExtension", ...
        "Selected %s file must end in %s; got: %s", ...
        role, strjoin(extensions, " / "), selectedPath);
end
if ~isfile(selectedPath)
    error("radia:simulink:IHGeometryBrowseMissing", ...
        "Selected geometry file does not exist: %s", selectedPath);
end

selectedPath = string(java.io.File(char(selectedPath)).getCanonicalPath());
set_param(blockPath, parameterName, char(selectedPath));
end
