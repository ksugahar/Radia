function result = benchmark_beam_transfer_mex_matlab(outputPath, matlabRoot, casePath, options)
%BENCHMARK_BEAM_TRANSFER_MEX_MATLAB Measure the transfer map through MEX.
arguments
    outputPath (1,1) string
    matlabRoot (1,1) string
    casePath (1,1) string
    options.Repeats (1,1) double {mustBeInteger,mustBePositive} = 31
end

rawCase = fileread(casePath);
caseData = jsondecode(rawCase);
validateCase(caseData);
lengths = double(caseData.lengths_m(:));
A = permute(double(caseData.A_per_m), [2, 3, 1]);
F2 = zeros(6, 6, 6, numel(lengths));
F3 = zeros(6, 6, 6, 6, numel(lengths));
for index = 1:size(caseData.F2_entries, 1)
    item = caseData.F2_entries(index, :);
    F2(item(2)+1, item(3)+1, item(4)+1, item(1)+1) = item(5);
end
for index = 1:size(caseData.F3_entries, 1)
    item = caseData.F3_entries(index, :);
    F3(item(2)+1, item(3)+1, item(4)+1, item(5)+1, item(1)+1) = item(6);
end

addpath(matlabRoot);
setupInfo = radia.setup(Force=true, ConfigureSimulinkFileGeneration=false);
mexPath = string(setupInfo.mex_path);
if ~isfile(mexPath)
    error("radia:benchmark:MexMissing", "radia_mex not found at %s", mexPath);
end
repoRoot = string(fileparts(matlabRoot));
if ~strcmpi(string(fileparts(mexPath)), matlabRoot)
    error("radia:benchmark:MexSource", ...
        "radia_mex resolved outside the selected MATLAB source tree: %s", mexPath);
end
requireCleanSource(repoRoot);
buildManifestPath = mexPath + ".build.json";
buildProvenance = readBuildProvenance(buildManifestPath, mexPath, repoRoot);
sourceCommit = string(buildProvenance.source_commit);
mexFile = dir(mexPath);
api = radia.apiInfo();
sourceVersion = radiaSourceVersion(repoRoot);
run = @() radia.beam.propagateVariationalMap( ...
    lengths, A, F2PerM=F2, F3PerM=F3, Names=string(caseData.names), ...
    MaximumOrder=caseData.maximum_order, ...
    MaximumStepM=caseData.maximum_step_m);

timer = tic;
first = run();
firstSeconds = toc(timer);
samples = zeros(options.Repeats, 1);
for index = 1:options.Repeats
    timer = tic;
    sample = run(); %#ok<NASGU>
    samples(index) = toc(timer);
end
if ~isfinite(firstSeconds) || firstSeconds < 0 || ...
        any(~isfinite(samples), "all") || any(samples <= 0, "all")
    error("radia:benchmark:Timing", "Benchmark timings must be finite and positive");
end
measuredObservables = observables(first);
if any(~isfinite(cell2mat(struct2cell(measuredObservables))), "all")
    error("radia:benchmark:Observables", "Benchmark observables must be finite");
end

result = struct( ...
    "schema", "radia.beam-transfer-backend-benchmark/v1", ...
    "executed_at_utc", datetime("now", TimeZone="UTC", ...
        Format="yyyy-MM-dd'T'HH:mm:ss.SSSXXX"), ...
    "backend", "matlab-mex", ...
    "case_id", string(caseData.case_id), ...
    "case_sha256", sha256File(casePath), ...
    "radia_version", sourceVersion, ...
    "radia_legacy_core_version", radia.UtiVer(), ...
    "radia_git_head", sourceCommit, ...
    "matlab_version", string(version), ...
    "host", string(getenv("COMPUTERNAME")), ...
    "binary", struct( ...
        "path", "matlab/" + string(mexFile.name), ...
        "bytes", mexFile.bytes, ...
        "sha256", sha256File(mexPath), ...
        "source_commit", sourceCommit, ...
        "build_manifest", "matlab/" + string(mexFile.name) + ".build.json", ...
        "modified", datetime(mexFile.datenum, ConvertFrom="datenum", ...
            TimeZone="local", Format="yyyy-MM-dd'T'HH:mm:ssXXX"), ...
        "api_version", api.api_version), ...
    "repeats", options.Repeats, ...
    "first_s", firstSeconds, ...
    "median_s", median(samples), ...
    "min_s", min(samples), ...
    "observables", measuredObservables);

outputDirectory = fileparts(outputPath);
if strlength(outputDirectory) > 0 && ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
file = fopen(outputPath, "w");
if file < 0
    error("radia:benchmark:Output", "Could not open %s", outputPath);
end
cleanup = onCleanup(@() fclose(file));
fwrite(file, jsonencode(result, PrettyPrint=true), "char");
fprintf("%s\n", jsonencode(result, PrettyPrint=true));
end

function validateCase(value)
expectedSchema = "radia.beam-transfer-benchmark-case/v1";
if ~isfield(value, "schema") || string(value.schema) ~= expectedSchema
    error("radia:benchmark:CaseSchema", "Case schema must be %s", expectedSchema);
end
lengths = double(value.lengths_m(:));
if isempty(lengths) || any(~isfinite(lengths)) || any(lengths <= 0)
    error("radia:benchmark:CaseLengths", "lengths_m must be finite and positive");
end
count = numel(lengths);
if ~isequal(size(value.A_per_m), [count, 6, 6]) || any(~isfinite(value.A_per_m), "all")
    error("radia:benchmark:CaseA", "A_per_m must be finite with shape [segment,6,6]");
end
if numel(value.names) ~= count
    error("radia:benchmark:CaseNames", "names must match lengths_m");
end
if ~ismember(value.maximum_order, [1, 2, 3])
    error("radia:benchmark:CaseOrder", "maximum_order must be 1, 2, or 3");
end
if ~isscalar(value.maximum_step_m) || ~isfinite(value.maximum_step_m) || value.maximum_step_m <= 0
    error("radia:benchmark:CaseStep", "maximum_step_m must be finite and positive");
end
validateSparseEntries(value.F2_entries, count, 5, "F2_entries");
validateSparseEntries(value.F3_entries, count, 6, "F3_entries");
end

function validateSparseEntries(entries, count, width, name)
if isempty(entries)
    return
end
if size(entries, 2) ~= width || any(~isfinite(entries), "all")
    error("radia:benchmark:CaseEntries", "%s has an invalid shape or value", name);
end
indices = entries(:, 1:end-1);
limits = [count, repmat(6, 1, width-2)];
if any(indices ~= fix(indices), "all") || any(indices < 0, "all") || ...
        any(indices >= limits, "all")
    error("radia:benchmark:CaseEntries", "%s contains an out-of-range index", name);
end
end

function result = observables(map)
result = struct( ...
    "R_fro", norm(map.R, "fro"), ...
    "T_fro", norm(map.T(:)), ...
    "U_fro", norm(map.U(:)), ...
    "R_11", map.R(1,1), ...
    "T_211", map.T(2,1,1), ...
    "U_3111", map.U(3,1,1,1), ...
    "U_4111", map.U(4,1,1,1), ...
    "R_composition_error", map.diagnostics.R_composition_error, ...
    "T_reconstruction_error", map.diagnostics.T_reconstruction_error, ...
    "U_reconstruction_error", map.diagnostics.U_reconstruction_error);
end

function value = gitHead(repoRoot)
[status, output] = system('git -C "' + string(repoRoot) + '" rev-parse HEAD');
if status ~= 0
    error("radia:benchmark:GitHead", "Could not resolve source commit for %s", repoRoot);
end
value = strtrim(string(output));
if strlength(value) ~= 40
    error("radia:benchmark:GitHead", "Expected a full source commit, got %s", value);
end
end

function requireCleanSource(repoRoot)
[workingStatus, ~] = system('git -C "' + string(repoRoot) + '" diff --quiet HEAD --');
[stagedStatus, ~] = system('git -C "' + string(repoRoot) + '" diff --cached --quiet HEAD --');
if workingStatus ~= 0 || stagedStatus ~= 0
    error("radia:benchmark:DirtySource", ...
        "Beam backend validation requires a clean source checkout: %s", repoRoot);
end
end

function manifest = readBuildProvenance(manifestPath, binaryPath, repoRoot)
if ~isfile(manifestPath)
    error("radia:benchmark:BuildProvenance", ...
        "Native build provenance manifest is missing: %s", manifestPath);
end
manifest = jsondecode(fileread(manifestPath));
expectedSchema = "radia.native-build-provenance.v1";
if ~isfield(manifest, "schema") || string(manifest.schema) ~= expectedSchema
    error("radia:benchmark:BuildProvenance", ...
        "Native provenance schema must be %s", expectedSchema);
end
if ~isfield(manifest, "source_dirty") || ~isequal(manifest.source_dirty, false)
    error("radia:benchmark:BuildProvenance", ...
        "Native provenance must come from a clean source build");
end
if ~isfield(manifest, "source_change_fingerprint_sha256") || ...
        strlength(string(manifest.source_change_fingerprint_sha256)) ~= 64
    error("radia:benchmark:BuildProvenance", ...
        "Native provenance source change fingerprint is invalid");
end
sourceCommit = string(manifest.source_commit);
if sourceCommit ~= gitHead(repoRoot)
    error("radia:benchmark:BuildProvenance", ...
        "Native provenance source commit does not match the selected checkout");
end
binary = dir(binaryPath);
if string(manifest.binary_name) ~= string(binary.name) || ...
        manifest.binary_bytes ~= binary.bytes || ...
        string(manifest.binary_sha256) ~= sha256File(binaryPath)
    error("radia:benchmark:BuildProvenance", ...
        "Native provenance does not match the loaded binary");
end
end

function value = radiaSourceVersion(repoRoot)
projectFile = fullfile(repoRoot, "pyproject.toml");
if ~isfile(projectFile)
    error("radia:benchmark:SourceVersion", "Missing %s", projectFile);
end
match = regexp(fileread(projectFile), ...
    '(?m)^version\s*=\s*"([^"]+)"\s*$', "tokens", "once");
if isempty(match)
    error("radia:benchmark:SourceVersion", ...
        "Could not read Radia version from %s", projectFile);
end
value = string(match{1});
end

function value = sha256File(path)
file = fopen(path, "rb");
if file < 0
    error("radia:benchmark:HashInput", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, Inf, "*uint8");
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(bytes);
hashBytes = typecast(digest.digest(), "uint8");
value = lower(string(reshape(dec2hex(hashBytes, 2).', 1, [])));
end
