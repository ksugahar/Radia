function inside = installed_wheel_path_inside(file, directory)
%INSTALLED_WHEEL_PATH_INSIDE True when file lies inside directory.
%
% Both paths are compared after separator, trailing-separator and (on
% Windows) case normalization, with a directory-separator boundary, so a
% sibling directory that shares the prefix is outside. Relative segments
% fail loudly instead of being resolved.

arguments
    file (1,1) string
    directory (1,1) string
end
file = normalize(file);
directory = normalize(directory);
inside = file ~= "" && directory ~= "" && startsWith(file, directory + "/");
end

function value = normalize(value)
value = replace(strtrim(value), "\", "/");
while strlength(value) > 1 && endsWith(value, "/")
    value = extractBefore(value, strlength(value));
end
segments = split(value, "/");
if any(segments == "." | segments == "..")
    error("radia:optuna:InstalledWheelPath", ...
        "Relative path segments are not accepted: %s", value);
end
if ispc
    value = lower(value);
end
end
