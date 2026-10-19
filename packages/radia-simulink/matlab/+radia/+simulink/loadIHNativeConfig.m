function config = loadIHNativeConfig(configFile)
%LOADIHNATIVECONFIG Load an IH native S-Function configuration.

arguments
    configFile {mustBeTextScalar} = ""
end

configFile = string(configFile);
if strlength(configFile) == 0
    if evalin("base", "exist('radia_ih_config','var')") ~= 1
        error("radia:simulink:IHConfigMissing", ...
            "Define radia_ih_config or select a native IH MAT/JSON configuration.");
    end
    config = evalin("base", "radia_ih_config");
else
    if ~isfile(configFile)
        error("radia:simulink:IHConfigFile", ...
            "Native IH configuration does not exist: %s", configFile);
    end
    [~, ~, extension] = fileparts(configFile);
    switch lower(extension)
        case ".mat"
            payload = load(configFile);
        case ".json"
            payload = jsondecode(fileread(configFile));
        otherwise
            error("radia:simulink:IHConfigFileType", ...
                "Native IH configuration must be a MAT or JSON file.");
    end
    if isfield(payload, "config")
        config = payload.config;
    elseif isfield(payload, "radia_ih_config")
        config = payload.radia_ih_config;
    else
        config = payload;
    end
end
if ~isstruct(config) || ~isscalar(config)
    error("radia:simulink:IHConfigStruct", ...
        "Native IH configuration must be one scalar struct.");
end
end
