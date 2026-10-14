function canonical = canonicalInputPath(path)
%CANONICALINPUTPATH Canonical absolute spelling of an existing input file.
% Separators, "." and ".." segments and letter case follow the file system,
% matching Python's Path.resolve() for the provenance record.
arguments
    path (1,1) string
end
file = java.io.File(char(path));
if ~file.isAbsolute() || ~isfile(path)
    error("radia:simulink:IHInputPath", "Input requires an existing absolute file path: %s", path);
end
canonical = string(file.getCanonicalPath());
end
