function status = updateIHGeometry(modelName, options)
%UPDATEIHGEOMETRY Rebuild IH operators when the geometry files changed.
%   status = updateIHGeometry(modelName) drives the "Geometry Update"
%   block (Tag RadiaIHGeometryUpdate, see addIHGeometryUpdateBlock):
%   it reads the block's workpiece/coil file paths, fingerprints their
%   CONTENT (SHA-256) together with the assemble command, and compares
%   against the sidecar written next to the configuration file.  When
%   anything changed, the assemble command is executed (this is the
%   explicit-update boundary where expensive mesh/basis/matrix
%   construction is allowed), the refreshed configuration is loaded
%   through radia.simulink.configureIHNativeModel, and the sidecar is
%   rewritten with an incremented revision.  When nothing changed the
%   configuration is only (re)loaded -- no command runs.
%
%   The model InitFcn installed by addIHGeometryUpdateBlock calls this
%   on every diagram update / simulation start, so re-pointing or
%   editing the .vol/.step files is picked up automatically at the next
%   update; the mask's "Rebuild now" button calls it with Force=true.
%
%   status fields: engaged, rebuilt, revision, reason, notes,
%   config_file, files.
%
%   Fail-loud contract (no fallbacks): a missing geometry file, a
%   failing assemble command, a command that does not produce the
%   configuration file, and a stale state with auto-rebuild disabled
%   are all immediate errors.

arguments
    modelName (1,1) string
    options.Force (1,1) logical = false
end

status = struct("engaged", false, "rebuilt", false, "revision", 0, ...
    "reason", "", "notes", strings(0, 1), "config_file", "", ...
    "files", strings(0, 1));

blocks = find_system(modelName, "LookUnderMasks", "all", ...
    "FollowLinks", "on", "Tag", "RadiaIHGeometryUpdate");
if isempty(blocks)
    status.reason = "no-geometry-update-block";
    return
end
if numel(blocks) > 1
    error("radia:simulink:IHGeometryUpdateDuplicate", ...
        "Model %s contains %d Geometry Update blocks; keep exactly one.", ...
        modelName, numel(blocks));
end
block = blocks{1};

wpValue = strtrim(string(get_param(block, "wp_vol")));
coilValue = strtrim(string(get_param(block, "coil_file")));
command = strtrim(string(get_param(block, "assemble_command")));
configFile = strtrim(string(get_param(block, "config_file")));
autoRebuild = strcmpi(string(get_param(block, "auto_rebuild")), "on");

if strlength(wpValue) == 0 && strlength(coilValue) == 0 && ...
        strlength(command) == 0 && strlength(configFile) == 0
    % Block placed but not configured yet: nothing to watch.
    status.reason = "unconfigured";
    return
end
if strlength(wpValue) == 0 || strlength(coilValue) == 0
    error("radia:simulink:IHGeometryUpdateInputs", ...
        "Geometry Update needs both a workpiece file and a coil file.");
end
if strlength(configFile) == 0
    error("radia:simulink:IHGeometryUpdateConfig", ...
        "Geometry Update needs the configuration file the assemble " + ...
        "command writes (config_file).");
end

[wpPath, coilRole, coilPath, notes] = classifyGeometryPair(wpValue, coilValue);
status.engaged = true;
status.notes = notes;
status.config_file = char(configFile);
status.files = [wpPath; coilPath];
for k = 1:numel(notes)
    warning("radia:simulink:IHGeometryRolesReassigned", "%s", notes(k));
end

fingerprint = struct( ...
    "schema", "radia.ih.simulink.geometry_fingerprint.v1", ...
    "command", char(command), ...
    "coil_role", char(coilRole), ...
    "files", radia.simulink.fileFingerprint([wpPath; coilPath]));

sidecarPath = configFile + ".fingerprint.json";
previousRevision = 0;
fresh = false;
if ~options.Force && isfile(sidecarPath) && isfile(configFile)
    stored = jsondecode(fileread(sidecarPath));
    previousRevision = storedRevision(stored);
    fresh = sameFingerprint(stored, fingerprint);
elseif isfile(sidecarPath)
    stored = jsondecode(fileread(sidecarPath));
    previousRevision = storedRevision(stored);
end

if fresh
    radia.simulink.configureIHNativeModel(modelName, configFile);
    status.revision = previousRevision;
    status.reason = "up-to-date";
    return
end

if ~autoRebuild
    error("radia:simulink:IHGeometryUpdateStale", ...
        "Geometry inputs changed but auto rebuild is off.  Press " + ...
        "'Rebuild now' on the Geometry Update block or enable " + ...
        "auto_rebuild.\n  workpiece: %s\n  coil (%s): %s", ...
        wpPath, coilRole, coilPath);
end
if strlength(command) == 0
    error("radia:simulink:IHGeometryUpdateCommand", ...
        "Geometry inputs changed but no assemble command is set on " + ...
        "the Geometry Update block.");
end

[exitCode, output] = system(char(command));
if exitCode ~= 0
    error("radia:simulink:IHGeometryUpdateAssemble", ...
        "Assemble command failed (exit %d):\n  %s\n--- output tail ---\n%s", ...
        exitCode, command, tailOf(output, 2000));
end
if ~isfile(configFile)
    error("radia:simulink:IHGeometryUpdateArtifact", ...
        "Assemble command exited 0 but did not write the " + ...
        "configuration file: %s", configFile);
end
radia.simulink.configureIHNativeModel(modelName, configFile);

record = fingerprint;
record.revision = previousRevision + 1;
record.generated_by = "radia.simulink.updateIHGeometry";
writeJson(sidecarPath, record);

status.rebuilt = true;
status.revision = record.revision;
status.reason = "rebuilt";
end

function [wpPath, coilRole, coilPath, notes] = classifyGeometryPair(a, b)
%CLASSIFYGEOMETRYPAIR Assign the two mask values to roles by extension.
%   The workpiece slot needs a mesh (.vol/.vol.gz); the coil slot takes
%   a STEP (.step/.stp -> PEEC filaments) or a mesh (.vol -> BEM-A).
%   A crossed pair (STEP typed into the workpiece box) is repaired
%   deterministically and reported; two STEP files cannot form a valid
%   pair and error immediately.
notes = strings(0, 1);
volExtensions = [".vol", ".vol.gz"];
stepExtensions = [".step", ".stp"];
aIsVol = fitsExtension(a, volExtensions);
bIsVol = fitsExtension(b, volExtensions);
aIsStep = fitsExtension(a, stepExtensions);
bIsStep = fitsExtension(b, stepExtensions);
if ~((aIsVol || aIsStep) && (bIsVol || bIsStep))
    bad = a; if aIsVol || aIsStep, bad = b; end
    error("radia:simulink:IHGeometryUpdateExtension", ...
        "Geometry Update accepts .vol/.vol.gz meshes and .step/.stp " + ...
        "CAD only; got: %s", bad);
end
if aIsVol && bIsStep
    wpPath = a; coilRole = "peec_step"; coilPath = b;
elseif aIsStep && bIsVol
    wpPath = b; coilRole = "peec_step"; coilPath = a;
    notes(end + 1, 1) = sprintf("geometry input reassigned by " + ...
        "extension: wp_vol <- '%s' (was '%s')", b, a);
    notes(end + 1, 1) = sprintf("geometry input reassigned by " + ...
        "extension: peec_step <- '%s' (was '%s')", a, b);
elseif aIsVol && bIsVol
    wpPath = a; coilRole = "coil_vol"; coilPath = b;
else
    error("radia:simulink:IHGeometryUpdateTwoSteps", ...
        "Both geometry inputs are STEP files; the workpiece must be " + ...
        "a .vol mesh:\n  %s\n  %s", a, b);
end
end

function tf = fitsExtension(value, extensions)
tf = any(endsWith(lower(string(value)), lower(string(extensions))));
end

function revision = storedRevision(stored)
revision = 0;
if isfield(stored, "revision") && isscalar(stored.revision) && ...
        isfinite(stored.revision)
    revision = double(stored.revision);
end
end

function tf = sameFingerprint(stored, current)
tf = false;
if ~isfield(stored, "command") || ~isfield(stored, "files") || ...
        ~isfield(stored, "coil_role")
    return
end
if ~strcmp(char(string(stored.command)), current.command) || ...
        ~strcmp(char(string(stored.coil_role)), current.coil_role)
    return
end
storedFiles = stored.files;
if numel(storedFiles) ~= numel(current.files)
    return
end
for index = 1:numel(current.files)
    storedEntry = storedFiles(index);
    if ~isfield(storedEntry, "sha256") || ...
            ~strcmp(char(string(storedEntry.sha256)), ...
                current.files(index).sha256)
        return
    end
end
tf = true;
end

function writeJson(path, record)
identifier = fopen(path, "w");
if identifier < 0
    error("radia:simulink:IHGeometryUpdateSidecar", ...
        "Cannot write fingerprint sidecar: %s", path);
end
closer = onCleanup(@() fclose(identifier));
fwrite(identifier, jsonencode(record, PrettyPrint=true), "char");
end

function text = tailOf(text, limit)
text = char(text);
if numel(text) > limit
    text = text(end - limit + 1:end);
end
end
