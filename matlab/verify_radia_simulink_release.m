function report = verify_radia_simulink_release()
%VERIFY_RADIA_SIMULINK_RELEASE Verify an extracted full Radia library package.
matlabRoot = fileparts(mfilename("fullpath"));
pathWasPresent = isPathEntry(matlabRoot);
if ~pathWasPresent
    addpath(matlabRoot, "-begin");
end
cleanupPath = onCleanup(@() removePathIfAdded(matlabRoot, pathWasPresent));

required = ["radia_simulink_library.slx", "radia_ih.slx", ...
    "radia_streamfunction_optimization.slx", "install_radia_simulink.m", ...
    "radia_mex." + mexext, "radia_ih_eddy_sfun.m", ...
    "radia_ih_thermal_sfun.m"];
missing = required(~isfile(fullfile(matlabRoot, required)));
if ~isempty(missing)
    error("radia:simulink:ReleaseMissing", ...
        "Radia Simulink release is missing: %s", strjoin(missing, ", "));
end

setupInfo = radia.setup(RequireMex=true,Force=true);
if ~setupInfo.mex_available
    error("radia:simulink:ReleaseMex", ...
        "The full library package did not load radia_mex.");
end
api = radia.apiInfo();

libraryPath = fullfile(matlabRoot, "radia_simulink_library.slx");
load_system(libraryPath);
cleanupLibrary = onCleanup(@() closeIfLoaded("radia_simulink_library"));
if string(get_param("radia_simulink_library", "BlockDiagramType")) ~= "library"
    error("radia:simulink:ReleaseLibrary", ...
        "radia_simulink_library.slx is not a Simulink library.");
end
motorPath = ...
    "radia_simulink_library/Reduced Models/Motor Angle Family";
if getSimulinkBlockHandle(motorPath) <= 0 || ...
        string(get_param(motorPath, "FunctionName")) ~= ...
        "radia_motor_angle_family_mex_sfunction"
    error("radia:simulink:ReleaseMotorBlock", ...
        "The native Motor Angle Family block is absent or miswired.");
end
motorContract = get_param(motorPath, "UserData");
if string(motorContract.backend) ~= ...
        "native-mex-periodic-interpolation" || motorContract.python_per_step
    error("radia:simulink:ReleaseMotorBackend", ...
        "The Motor Angle Family block backend contract is invalid.");
end

grid = [0; pi];
period = 2*pi;
A = reshape([1.0, 0.5], 1, 1, 2);
B = reshape([1.0, 3.0], 1, 1, 2);
C = reshape([1.0, 3.0], 1, 1, 2);
D = zeros(1, 1, 2);
Q = reshape([2.0, 4.0], 1, 1, 2);
R = reshape([1.0, 2.0], 1, 1, 2);
S = reshape([0.0, 2.0], 1, 1, 2);
nativeHandle = radia.internal.callMex( ...
    "simulink.state_space.create", grid, period, A, B, C, D, ...
    Q, R, S, 2.0);
cleanupHandle = onCleanup(@() destroyNativeStateSpace(nativeHandle));
first = radia.internal.callMex( ...
    "simulink.state_space.output", nativeHandle, pi/2, 3.0);
radia.internal.callMex( ...
    "simulink.state_space.update", nativeHandle, pi/2, 3.0);
second = radia.internal.callMex( ...
    "simulink.state_space.output", nativeHandle, pi/2, 3.0);
snapshot = radia.internal.callMex( ...
    "simulink.state_space.snapshot", nativeHandle);
if norm(first - [4.0; 19.5], inf) > 1e-12 || ...
        norm(second - [15.0; 122.625], inf) > 1e-12 || ...
        snapshot.step_count ~= 1
    error("radia:simulink:ReleaseMotorNumerics", ...
        "The extracted native motor lifecycle smoke failed.");
end
clear cleanupHandle

ih = verify_radia_ih_release();
report = struct( ...
    "passed", true, ...
    "matlab_release", string(version("-release")), ...
    "mex_api_version", api.api_version, ...
    "library", "radia_simulink_library", ...
    "motor_backend", string(motorContract.backend), ...
    "motor_first_output", first, ...
    "motor_second_output", second, ...
    "ih_backend", string(ih.backend), ...
    "python_per_step", false);
fprintf("RADIA_SIMULINK_RELEASE_OK matlab=%s motor=%s ih=%s\n", ...
    report.matlab_release, report.motor_backend, report.ih_backend);
clear cleanupLibrary cleanupPath
end

function destroyNativeStateSpace(handle)
try
    radia.internal.callMex("simulink.state_space.destroy", handle);
catch
end
end

function closeIfLoaded(name)
if bdIsLoaded(name)
    close_system(name, 0);
end
end

function present = isPathEntry(folder)
present = any(strcmpi(split(string(path), pathsep), string(folder)));
end

function removePathIfAdded(folder, pathWasPresent)
if ~pathWasPresent && isPathEntry(folder)
    rmpath(folder);
end
end
