function onMagLevBlockFamilyChanged(blockPath)
%ONMAGLEVBLOCKFAMILYCHANGED Apply a MagLev plant family expression.

arguments
    blockPath (1,1) string
end

expression = string(get_param(blockPath, "family"));
if strlength(strtrim(expression)) == 0
    error("radia:simulink:MagLevFamily", ...
        "Common-basis CLN family expression must not be empty.");
end
% Force Simulink to materialize the mask workspace before changing a child
% S-Function parameter. Without this call, set_param evaluates the new
% expression against an uninitialized mask workspace and rejects even valid
% package functions and model-workspace variables.
mask = Simulink.Mask.get(blockPath);
if isempty(mask)
    error("radia:simulink:MagLevMask", ...
        "The MagLev plant must remain a masked subsystem.");
end
mask.getWorkspaceVariables;
set_param(blockPath + "/Moving HCurl CLN", ...
    "Parameters", char(expression));
end
