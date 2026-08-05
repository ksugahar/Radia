function status = updateIHGeometry(modelName, options)
%UPDATEIHGEOMETRY Rebuild IH operators when the geometry files changed.
%   status = updateIHGeometry(modelName) drives the "Geometry Update"
%   block (Tag RadiaIHGeometryUpdate, see addIHGeometryUpdateBlock):
%   it reads the block's workpiece/coil file paths, fingerprints their
%   CONTENT (SHA-256) together with the configured assembler, and
%   compares against the sidecar written next to the configuration
%   file.  When anything changed, the assembler runs (this is the
%   explicit-update boundary where expensive mesh/basis/matrix
%   construction is allowed), the refreshed configuration is loaded
%   through radia.simulink.configureIHNativeModel, and the sidecar is
%   rewritten with an incremented revision.  When nothing changed the
%   configuration is only (re)loaded -- no assembler runs.
%
%   With both custom assembler fields empty, the built-in
%   assembleIHOperatorsFromGeometry function runs Radia's shape-to-
%   operator CLI.  Otherwise the assembler is EITHER assemble_fcn --
%   called in-process as fcn(wpVol, coilFile, configFile) -- OR
%   assemble_command, a shell command. Setting both is an error.
%
%   The model InitFcn installed by addIHGeometryUpdateBlock calls this
%   on every diagram update / simulation start, so re-pointing or
%   editing the .vol/.step files is picked up automatically at the next
%   update; the mask's "Rebuild now" button calls it with Force=true.
%
%   status fields: engaged, rebuilt, revision, reason, notes,
%   config_file, files.  reason is one of:
%     "no-geometry-update-block" / "unconfigured"  (nothing to do)
%     "up-to-date"  inputs unchanged AND the workspace provably holds
%                   this configuration already (artifact hash +
%                   revision marker) -- nothing ran
%     "reloaded"    inputs unchanged but the configuration file content
%                   changed (or the workspace lost it) -- reloaded
%                   without running the assemble command
%     "rebuilt"     inputs changed -- assemble command ran
%
%   Geometry paths must be ABSOLUTE. A blank config_file is derived
%   beside the workpiece from the two geometry names.
%
%   Fail-loud contract (no fallbacks): a relative path, a missing
%   geometry file, a failing assemble command, a command that does not
%   produce the configuration file, and a stale state with auto-rebuild
%   disabled (unless Force=true was explicitly requested) are all
%   immediate errors.

arguments
    modelName (1,1) string
    options.Force (1,1) logical = false
end

status = struct("engaged", false, "rebuilt", false, "revision", 0, ...
    "reason", "", "notes", strings(0, 1), "config_file", "", ...
    "files", strings(0, 1), "assembler", "");

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
assembleFcn = strtrim(string(get_param(block, "assemble_fcn")));
command = strtrim(string(get_param(block, "assemble_command")));
configFile = strtrim(string(get_param(block, "config_file")));
autoRebuild = strcmpi(string(get_param(block, "auto_rebuild")), "on");

if strlength(wpValue) == 0 && strlength(coilValue) == 0 && ...
        strlength(configFile) == 0
    % Block placed but not configured yet: nothing to watch.
    status.reason = "unconfigured";
    return
end
if strlength(wpValue) == 0 || strlength(coilValue) == 0
    error("radia:simulink:IHGeometryUpdateInputs", ...
        "Geometry Update needs both a workpiece file and a coil file.");
end
% The InitFcn hook fires from whatever working directory MATLAB happens
% to be in, so relative paths would resolve differently between
% sessions (wrong file hashed, sidecar scattered).  Require absolute.
pathValues = [wpValue, coilValue];
for pathIndex = 1:numel(pathValues)
    if ~java.io.File(char(pathValues(pathIndex))).isAbsolute()
        error("radia:simulink:IHGeometryUpdateRelativePath", ...
            "Geometry Update requires ABSOLUTE paths (the update hook " + ...
            "runs from an arbitrary working directory); got: %s", ...
            pathValues(pathIndex));
    end
end

[wpPath, coilRole, coilPath, notes] = classifyGeometryPair(wpValue, coilValue);
if strlength(configFile) == 0
    configFile = derivedConfigFile(wpPath, coilPath);
    set_param(block, "config_file", char(configFile));
elseif ~java.io.File(char(configFile)).isAbsolute()
    error("radia:simulink:IHGeometryUpdateRelativePath", ...
        "Geometry Update requires an ABSOLUTE config_file path; got: %s", ...
        configFile);
end
assemblyOptions = readAssemblyOptions(block);
hasFcn = strlength(assembleFcn) > 0;
hasCommand = strlength(command) > 0;
if hasFcn && hasCommand
    error("radia:simulink:IHGeometryUpdateAmbiguousAssemble", ...
        "Set either assemble_fcn or assemble_command, not both.");
end
builtInFcn = "radia.simulink.assembleIHOperatorsFromGeometry";
if ~hasFcn && ~hasCommand
    assembleFcn = builtInFcn;
    hasFcn = true;
end
status.engaged = true;
status.notes = notes;
status.config_file = char(configFile);
status.files = [wpPath; coilPath];
for k = 1:numel(notes)
    warning("radia:simulink:IHGeometryRolesReassigned", "%s", notes(k));
end

fingerprint = struct( ...
    "schema", "radia.ih.simulink.geometry_fingerprint.v2", ...
    "assemble_fcn", char(assembleFcn), ...
    "command", char(command), ...
    "coil_role", char(coilRole), ...
    "assembly_options", assemblyOptions, ...
    "files", radia.simulink.fileFingerprint([wpPath; coilPath]));

sidecarPath = configFile + ".fingerprint.json";
stored = readSidecar(sidecarPath);
previousRevision = 0;
fresh = false;
if ~isempty(stored)
    previousRevision = storedRevision(stored);
    if ~options.Force && isfile(configFile)
        fresh = sameFingerprint(stored, fingerprint);
    end
end

if fresh
    % Inputs unchanged -> never re-run the assemble command.  Reloading
    % the configuration is still skipped only when BOTH the artifact
    % hash and the model-workspace revision marker prove the workspace
    % already holds exactly this configuration (a real config carries
    % dense operator matrices -- reloading it on every diagram update
    % costs seconds and doubles peak memory for nothing).
    workspace = get_param(modelName, "ModelWorkspace");
    artifact = radia.simulink.fileFingerprint(configFile);
    if canSkipReload(workspace, stored, artifact, previousRevision)
        status.revision = previousRevision;
        status.reason = "up-to-date";
        return
    end
    radia.simulink.configureIHNativeModel(modelName, configFile);
    workspace.assignin("radia_ih_geometry_revision", previousRevision);
    record = stored;
    record.artifact = artifact;
    writeJson(sidecarPath, record);
    status.revision = previousRevision;
    status.reason = "reloaded";
    return
end

if ~autoRebuild && ~options.Force
    error("radia:simulink:IHGeometryUpdateStale", ...
        "Geometry inputs changed but auto rebuild is off.  Press " + ...
        "'Rebuild now' on the Geometry Update block or enable " + ...
        "auto_rebuild.\n  workpiece: %s\n  coil (%s): %s", ...
        wpPath, coilRole, coilPath);
end
if hasFcn
    % In-process MATLAB assembler (preferred): errors propagate with
    % their own stack, no shell quoting, and the function runs on the
    % MATLAB path like any other .m code.
    if assembleFcn == builtInFcn
        feval(char(assembleFcn), wpPath, coilPath, configFile, assemblyOptions);
        status.assembler = "radia-built-in-cli";
    else
        feval(char(assembleFcn), wpPath, coilPath, configFile);
        status.assembler = "matlab-function";
    end
else
    [exitCode, output] = system(char(command));
    if exitCode ~= 0
        error("radia:simulink:IHGeometryUpdateAssemble", ...
            "Assemble command failed (exit %d):\n  %s\n" + ...
            "--- output tail ---\n%s", ...
            exitCode, command, tailOf(output, 2000));
    end
    status.assembler = "shell-command";
end
if ~isfile(configFile)
    error("radia:simulink:IHGeometryUpdateArtifact", ...
        "The assembler finished but did not write the " + ...
        "configuration file: %s", configFile);
end
radia.simulink.configureIHNativeModel(modelName, configFile);

record = fingerprint;
record.artifact = radia.simulink.fileFingerprint(configFile);
record.revision = previousRevision + 1;
record.generated_by = "radia.simulink.updateIHGeometry";
writeJson(sidecarPath, record);
workspace = get_param(modelName, "ModelWorkspace");
workspace.assignin("radia_ih_geometry_revision", record.revision);

status.rebuilt = true;
status.revision = record.revision;
status.reason = "rebuilt";
end

function configFile = derivedConfigFile(wpPath, coilPath)
[folder, wpStem, wpExtension] = fileparts(wpPath);
if strcmpi(wpExtension, ".gz") && endsWith(lower(wpStem), ".vol")
    [~, wpStem] = fileparts(wpStem);
end
[~, coilStem, coilExtension] = fileparts(coilPath);
if strcmpi(coilExtension, ".gz") && endsWith(lower(coilStem), ".vol")
    [~, coilStem] = fileparts(coilStem);
end
wpStem = regexprep(string(wpStem), "[^A-Za-z0-9_.-]", "_");
coilStem = regexprep(string(coilStem), "[^A-Za-z0-9_.-]", "_");
configFile = string(fullfile(folder, ...
    wpStem + "_" + coilStem + "_ih_native.json"));
end

function values = readAssemblyOptions(block)
names = [ ...
    "python_executable", "frequency_hz", "coil_sigma", ...
    "workpiece_sigma", "workpiece_mu_r", "density", ...
    "heat_capacity", "thermal_conductivity", "convection", ...
    "initial_temperature_K", "sample_time_s", "workpiece_label", ...
    "coil_body_label", "coil_source_label", "coil_sink_label", ...
    "peec_n_peri", ...
    "peec_proximity", "coupling_mode", "workpiece_bem_backend"];
values = struct();
parameters = get_param(block, "ObjectParameters");
for index = 1:numel(names)
    name = names(index);
    if ~isfield(parameters, name)
        error("radia:simulink:IHGeometryUpdateLegacyBlock", ...
            "Geometry Update block lacks parameter '%s'. Rebuild the " + ...
            "model from radia.simulink.buildIHNativeModel.", name);
    end
    values.(name) = get_param(block, name);
end
end

function [wpPath, coilRole, coilPath, notes] = classifyGeometryPair(a, b)
%CLASSIFYGEOMETRYPAIR Assign the two mask values to roles by extension.
%   The workpiece slot needs a mesh (.vol/.vol.gz); the coil slot takes
%   a STEP (.step/.stp -> PEEC filaments) or a mesh (.vol -> BEM-A).
%   A crossed pair (STEP typed into the workpiece box) is repaired
%   deterministically and reported; two STEP files cannot form a valid
%   pair and error immediately.
notes = strings(0, 1);
extensions = radia.simulink.ihGeometryExtensions();
volExtensions = extensions.vol;
stepExtensions = extensions.step;
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

function stored = readSidecar(sidecarPath)
%READSIDECAR Sidecar record, or empty when missing or unreadable.
%   The sidecar is OUR cache metadata, not user input: a corrupted file
%   (e.g. a crash mid-write) must degrade to "no record" -- forcing a
%   rebuild from the real sources -- instead of blocking every diagram
%   update behind a cryptic jsondecode error.
stored = struct([]);
if ~isfile(sidecarPath)
    return
end
try
    stored = jsondecode(fileread(sidecarPath));
catch decodeError
    warning("radia:simulink:IHGeometryUpdateSidecarCorrupt", ...
        "Fingerprint sidecar is unreadable and will be rebuilt " + ...
        "(%s): %s", decodeError.message, sidecarPath);
    stored = struct([]);
end
if ~isstruct(stored) || ~isscalar(stored)
    stored = struct([]);
end
end

function tf = canSkipReload(workspace, stored, artifact, revision)
%CANSKIPRELOAD True when the model workspace provably already holds the
%   configuration the sidecar describes: the artifact hash matches the
%   config file on disk AND the workspace revision marker (assigned
%   together with every configure) matches the sidecar revision.  Any
%   doubt -> false -> reload (the safe direction).
tf = false;
if ~isfield(stored, "artifact")
    return
end
storedArtifact = stored.artifact;
if ~isfield(storedArtifact, "sha256") || ...
        ~strcmp(char(string(storedArtifact.sha256)), artifact.sha256)
    return
end
if ~workspace.hasVariable("radia_ih_config") || ...
        ~workspace.hasVariable("radia_ih_geometry_revision")
    return
end
marker = workspace.getVariable("radia_ih_geometry_revision");
if ~(isscalar(marker) && isnumeric(marker) && double(marker) == revision)
    return
end
tf = true;
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
% A sidecar written before the assemble_fcn field existed compares as
% different once, forcing one deterministic rebuild that upgrades it.
if ~isfield(stored, "schema") || ~isfield(stored, "command") || ...
        ~isfield(stored, "assemble_fcn") || ...
        ~isfield(stored, "assembly_options") || ...
        ~isfield(stored, "files") || ~isfield(stored, "coil_role")
    return
end
if ~strcmp(char(string(stored.schema)), current.schema) || ...
        ~strcmp(char(string(stored.command)), current.command) || ...
        ~strcmp(char(string(stored.assemble_fcn)), current.assemble_fcn) || ...
        ~strcmp(char(string(stored.coil_role)), current.coil_role) || ...
        ~strcmp(jsonencode(stored.assembly_options), ...
                jsonencode(current.assembly_options))
    return
end
storedFiles = stored.files;
if numel(storedFiles) ~= numel(current.files)
    return
end
for index = 1:numel(current.files)
    storedEntry = storedFiles(index);
    if ~isfield(storedEntry, "path") || ~isfield(storedEntry, "sha256") || ...
            ~strcmp(char(string(storedEntry.path)), ...
                current.files(index).path) || ...
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
