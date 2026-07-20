function configPath = writeApplicationConfig(application, settings, configPath, options)
%WRITEAPPLICATIONCONFIG Write a versioned DesignSpec configuration JSON.

arguments
    application (1,1) string {mustBeMember(application, ...
        ["em", "pcb", "motor", "streamfunction", "ih"])}
    settings (1,1) struct
    configPath (1,1) string
    options.PrimaryKey (1,1) string = ""
    options.WorkingDirectory (1,1) string = ""
end

parent = fileparts(configPath);
if strlength(parent) > 0 && ~isfolder(parent)
    mkdir(parent);
end

payload = struct( ...
    "schema", "radia.simulink.application_config.v1", ...
    "application", application, ...
    "settings", settings);
if strlength(options.PrimaryKey) > 0
    payload.primary_key = options.PrimaryKey;
end
if strlength(options.WorkingDirectory) > 0
    payload.working_directory = options.WorkingDirectory;
end

text = jsonencode(payload, PrettyPrint=true);
fileID = fopen(configPath, "w", "n", "UTF-8");
if fileID < 0
    error("radia:simulink:ConfigWriteFailed", ...
        "Could not open application config for writing: %s", configPath);
end
cleanup = onCleanup(@() fclose(fileID));
count = fprintf(fileID, "%s\n", text);
if count <= 0
    error("radia:simulink:ConfigWriteFailed", ...
        "Could not write application config: %s", configPath);
end
clear cleanup
configPath = string(configPath);
end
