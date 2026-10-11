function onIHConfigFileChanged(blockPath)
%ONIHCONFIGFILECHANGED Mask-callback entry point for config_file.
%   The .slx stores only "radia.simulink.onIHConfigFileChanged(gcb);"
%   -- all behavior lives in this named .m file (thin-callback policy:
%   inline mask-callback code is invisible to diffs and cannot be
%   tested).  Mask callbacks execute in the base workspace where mask
%   parameter variables are not defined, so the dialog value is read
%   back through get_param.

arguments
    blockPath (1,1) string
end

radia.simulink.configureIHNativeModel( ...
    string(bdroot(char(blockPath))), ...
    string(get_param(blockPath, "config_file")));
end
