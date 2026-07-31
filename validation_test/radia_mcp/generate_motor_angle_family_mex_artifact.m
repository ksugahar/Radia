% Generate the checked native Motor MEX/Simulink validation artifact.

scriptPath = string(mfilename("fullpath"));
root = string(fileparts(fileparts(fileparts(scriptPath))));
artifactDir = fullfile(root, "validation_test", "radia_mcp", "artifacts", ...
    "annular_motor_dual_lane_v1");
artifactPath = fullfile(artifactDir, "native_motor_angle_family.json");
matlabDir = fullfile(root, "matlab");
previousDirectory = string(pwd);
restoreDirectory = onCleanup(@() cd(previousDirectory));
cd(root);
addpath(matlabDir, "-begin");
clear radia.setup radia_mex

setupInfo = radia.setup(Force=true);
hostRole = string(getenv("RADIA_VALIDATION_HOST_ROLE"));
if strlength(hostRole) == 0
    hostRole = "developer-smoke";
end
if ~any(hostRole == ["compute","developer-smoke"])
    error("radia:validation:HostRole", ...
        "RADIA_VALIDATION_HOST_ROLE must be compute or developer-smoke.");
end
hostName = string(getenv("COMPUTERNAME"));
if strlength(hostName) == 0
    error("radia:validation:HostName", ...
        "COMPUTERNAME must identify the MATLAB validation host.");
end
testFiles = [ ...
    "tests/matlab/test_radia_mex.m", ...
    "tests/matlab/test_simulink_workflow.m"];
results = runtests(cellstr(fullfile(root, testFiles)));

testRows = repmat(struct( ...
    "name", "", "passed", false, "duration_s", 0.0), numel(results), 1);
for index = 1:numel(results)
    elapsed = results(index).Duration;
    if isduration(elapsed)
        elapsed = seconds(elapsed);
    end
    testRows(index) = struct( ...
        "name", string(results(index).Name), ...
        "passed", logical(results(index).Passed), ...
        "duration_s", double(elapsed));
end

excluded = string(setupInfo.excluded_openmp_runtime_dirs(:));
pathEntries = split(string(getenv("PATH")), pathsep);
excludedStillActive = false(size(excluded));
for index = 1:numel(excluded)
    excludedStillActive(index) = any(strcmpi(excluded(index), pathEntries));
end

sourcePath = fullfile(root, "src", "matlab", "radia_mex.cpp");
setupPath = fullfile(matlabDir, "+radia", "setup.m");
generatorPath = scriptPath + ".m";
mexPath = string(setupInfo.mex_path);
if ~startsWith(mexPath, root, IgnoreCase=ispc)
    error("radia:validation:MexLocation", ...
        "radia_mex must resolve inside the checkout under validation.");
end

artifact = struct( ...
    "radia_version", readProjectVersion(fullfile(root, "pyproject.toml")), ...
    "radia_core_api_version", string(radia.UtiVer()), ...
    "schema", "radia.validation.motor-angle-family-mex.v1", ...
    "executed_at_utc", utcTimestamp(), ...
    "execution_mode", "standalone_matlab_batch", ...
    "execution_policy", ...
        "MATLAB and Simulink verification with no external solver session.", ...
    "execution_environment", struct( ...
        "host_role", hostRole, ...
        "hostname", hostName, ...
        "platform", string(computer)), ...
    "matlab_version", string(version), ...
    "matlab_release", string(version("-release")), ...
    "simulink_version", toolboxVersion("Simulink"), ...
    "simulink_release", string(version("-release")), ...
    "optimization_toolbox_available", logical( ...
        ~isempty(ver("optim")) && ...
        license("test","Optimization_Toolbox") && ...
        exist("fmincon","file") == 2), ...
    "mex_api_version", double(radia.apiInfo().api_version), ...
    "mex_relative_path", relativePath(root, mexPath), ...
    "mex_sha256", sha256File(mexPath), ...
    "source_relative_path", relativePath(root, sourcePath), ...
    "source_sha256", sha256TextFile(sourcePath), ...
    "setup_relative_path", relativePath(root, setupPath), ...
    "setup_sha256", sha256TextFile(setupPath), ...
    "generator_relative_path", relativePath(root, generatorPath), ...
    "generator_sha256", sha256TextFile(generatorPath), ...
    "text_sha256_normalization", "newline-lf", ...
    "foreign_openmp_runtime_dirs_excluded_count", numel(excluded), ...
    "foreign_openmp_runtime_dirs_remaining_on_path_count", ...
        nnz(excludedStillActive), ...
    "test_files", testFiles, ...
    "test_count", numel(results), ...
    "passed_count", nnz([results.Passed]), ...
    "failed_count", nnz([results.Failed]), ...
    "incomplete_count", nnz([results.Incomplete]), ...
    "test_results", testRows, ...
    "validated_capabilities", [ ...
        "periodic_angle_family_native_interpolation", ...
        "quadratic_torque_output", ...
        "persistent_native_state", ...
        "split_output_update_lifecycle", ...
        "custom_sim_state_roundtrip", ...
        "simulink_s_function_compile", ...
        "foreign_openmp_runtime_isolation"], ...
    "status", statusFromResults(results));

if ~isfolder(artifactDir)
    mkdir(artifactDir);
end
writeUtf8Lf(artifactPath, ...
    string(jsonencode(artifact, PrettyPrint=true)) + string(newline));
fprintf("%s\n", jsonencode(struct( ...
    "status", artifact.status, ...
    "artifact", relativePath(root, artifactPath), ...
    "test_count", artifact.test_count, ...
    "passed_count", artifact.passed_count), PrettyPrint=true));
if artifact.status ~= "pass"
    error("radia:validation:MotorAngleFamily", ...
        "The native Motor MEX/Simulink validation suite did not pass.");
end
clear restoreDirectory

function value = readProjectVersion(path)
text = fileread(path);
token = regexp(text, '(?m)^version\s*=\s*"([^"]+)"', "tokens", "once");
if isempty(token)
    error("radia:validation:Version", "Could not read the project version.");
end
value = string(token{1});
end

function value = toolboxVersion(name)
details = ver(name);
if isempty(details)
    error("radia:validation:Toolbox", "%s is not installed.", name);
end
value = string(details(1).Version);
end

function value = utcTimestamp()
nowUtc = datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss.SSSXXX");
value = string(nowUtc);
end

function value = statusFromResults(results)
if ~isempty(results) && all([results.Passed]) && ...
        ~any([results.Failed]) && ~any([results.Incomplete])
    value = "pass";
else
    value = "fail";
end
end

function value = relativePath(root, path)
root = replace(string(root), "\", "/");
path = replace(string(path), "\", "/");
prefix = root + "/";
if ~startsWith(path, prefix, IgnoreCase=ispc)
    error("radia:validation:PublicPath", ...
        "Artifact paths must remain inside the checkout.");
end
value = extractAfter(path, strlength(prefix));
end

function value = sha256File(path)
file = fopen(path, "rb");
if file < 0
    error("radia:validation:HashRead", "Could not read %s.", path);
end
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, inf, "*uint8");
if usejava("jvm")
    digest = java.security.MessageDigest.getInstance("SHA-256");
    digest.update(typecast(bytes, "int8"));
    hashBytes = typecast(int8(digest.digest()), "uint8");
elseif ispc
    digest = System.Security.Cryptography.SHA256.Create();
    hashBytes = uint8(digest.ComputeHash(bytes));
else
    error("radia:validation:HashRuntime", ...
        "SHA-256 requires either the JVM or the Windows .NET runtime.");
end
value = lower(string(reshape(dec2hex(hashBytes, 2).', 1, [])));
end

function value = sha256TextFile(path)
file = fopen(path, "rb");
if file < 0
    error("radia:validation:HashRead", "Could not read %s.", path);
end
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, inf, "*uint8");
crlf = bytes(1:end-1) == 13 & bytes(2:end) == 10;
remove = [crlf; false];
bytes(remove) = [];
bytes(bytes == 13) = 10;
clear cleanup

if usejava("jvm")
    digest = java.security.MessageDigest.getInstance("SHA-256");
    digest.update(typecast(bytes, "int8"));
    hashBytes = typecast(int8(digest.digest()), "uint8");
elseif ispc
    digest = System.Security.Cryptography.SHA256.Create();
    hashBytes = uint8(digest.ComputeHash(bytes));
else
    error("radia:validation:HashRuntime", ...
        "SHA-256 requires either the JVM or the Windows .NET runtime.");
end
value = lower(string(reshape(dec2hex(hashBytes, 2).', 1, [])));
end

function writeUtf8Lf(path, text)
text = string(text);
if ~isscalar(text)
    error("radia:validation:ArtifactText", ...
        "The encoded artifact must be a scalar string.");
end
text = replace(string(text), sprintf("\r\n"), newline);
text = replace(text, sprintf("\r"), newline);
bytes = unicode2native(char(text), "UTF-8");
file = fopen(path, "wb");
if file < 0
    error("radia:validation:ArtifactWrite", "Could not write %s.", path);
end
cleanup = onCleanup(@() fclose(file));
fwrite(file, bytes, "uint8");
end
