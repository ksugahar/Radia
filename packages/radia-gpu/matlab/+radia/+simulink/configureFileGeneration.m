function info = configureFileGeneration(options)
%CONFIGUREFILEGENERATION Keep Simulink cache artifacts outside the repository.
%   INFO = radia.simulink.configureFileGeneration() configures Simulink's
%   cache and code-generation folders when they have not already been set.
%   Set RADIA_SIMULINK_FILEGEN_ROOT to override the default location.

arguments
    options.RootDirectory (1,1) string = ""
    options.Force (1,1) logical = false
    options.Verbose (1,1) logical = false
end

root = options.RootDirectory;
if strlength(root) == 0
    root = string(getenv("RADIA_SIMULINK_FILEGEN_ROOT"));
end
if strlength(root) == 0
    base = string(getenv("RADIA_TEMP_ROOT"));
    if strlength(base) == 0
        if ispc
            base = "C:\temp";
        else
            base = string(tempdir);
        end
    end
    root = fullfile(base, "radia", "simulink", ...
        "R" + string(version("-release")));
end

info = struct( ...
    "available", false, ...
    "changed", false, ...
    "root", root, ...
    "cache_folder", "", ...
    "codegen_folder", "", ...
    "reason", "simulink-unavailable");
if isempty(which("Simulink.fileGenControl")) || ...
        ~license("test", "Simulink")
    return
end

config = Simulink.fileGenControl("getConfig");
cacheFolder = string(config.CacheFolder);
codegenFolder = string(config.CodeGenFolder);
targetCache = fullfile(root, "cache");
targetCodegen = fullfile(root, "codegen");
sourceRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));

if options.Force
    cacheFolder = targetCache;
    codegenFolder = targetCodegen;
else
    if strlength(cacheFolder) == 0 || isInside(cacheFolder, sourceRoot)
        cacheFolder = targetCache;
    end
    if strlength(codegenFolder) == 0 || isInside(codegenFolder, sourceRoot)
        codegenFolder = targetCodegen;
    end
end

changed = string(config.CacheFolder) ~= cacheFolder || ...
    string(config.CodeGenFolder) ~= codegenFolder;
if changed
    Simulink.fileGenControl("set", ...
        CacheFolder=cacheFolder, ...
        CodeGenFolder=codegenFolder, ...
        createDir=true);
end

info = struct( ...
    "available", true, ...
    "changed", changed, ...
    "root", root, ...
    "cache_folder", cacheFolder, ...
    "codegen_folder", codegenFolder, ...
    "reason", "configured");
if options.Verbose
    fprintf("Radia Simulink generated files\n  Cache: %s\n  Codegen: %s\n", ...
        cacheFolder, codegenFolder);
end
end

function tf = isInside(candidate, parent)
candidate = string(java.io.File(candidate).getCanonicalPath());
parent = string(java.io.File(parent).getCanonicalPath());
if ispc
    candidate = lower(candidate);
    parent = lower(parent);
end
separator = string(filesep);
tf = candidate == parent || startsWith(candidate, parent + separator);
end
