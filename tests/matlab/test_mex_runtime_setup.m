function tests = test_mex_runtime_setup
%TEST_MEX_RUNTIME_SETUP Keep runtime checks without per-call Simulink mutation.
tests = functiontests(localfunctions);
end

function setupOnce(t)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root, "matlab"), "-begin");
t.TestData.Root = root;
radia.setup(Force=true);
end

function testCommandConversionAndRepeatedCalls(t)
expected = radia_mex('api.commands');
verifyEqual(t, radia.internal.callMex("api.commands"), expected);
verifyEqual(t, radia.internal.callMex('api.commands'), expected);
end

function testCallDoesNotConfigureSimulinkButExplicitSetupDoes(t)
original = Simulink.fileGenControl("getConfig");
cleanup = onCleanup(@() Simulink.fileGenControl("setConfig", config=original)); %#ok<NASGU>
Simulink.fileGenControl("set", CacheFolder=t.TestData.Root, CodeGenFolder=t.TestData.Root);
commands = radia.internal.callMex("api.commands");
verifyNotEmpty(t, commands);
current = Simulink.fileGenControl("getConfig");
verifyEqual(t, string(current.CacheFolder), string(t.TestData.Root));
verifyEqual(t, string(current.CodeGenFolder), string(t.TestData.Root));
info = radia.setup();
verifyTrue(t, info.simulink_file_generation.changed);
verifyNotEqual(t, info.simulink_file_generation.cache_folder, string(t.TestData.Root));
end

function testMissingMexIsCheckedAfterSuccessfulCall(t)
commands = radia.internal.callMex("api.commands");
verifyNotEmpty(t, commands);
folder = tempname("C:\temp");
mkdir(folder);
cleanup = onCleanup(@() removeShadow(folder)); %#ok<NASGU>
fid = fopen(fullfile(folder, "radia_mex.m"), "w");
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid, "function varargout=radia_mex(varargin)\nerror('test:unexpectedGateway','Missing-MEX check was bypassed');\nend\n");
clear fileCleanup
clear radia_mex
addpath(folder, "-begin");
verifyEqual(t, exist("radia_mex", "file"), 2);
verifyError(t, @() radia.internal.callMex("api.commands"), "radia:setup:MissingMex");
rmpath(folder);
clear radia_mex
verifyEqual(t, radia.internal.callMex("api.commands"), commands);
end

function testExplicitPythonChangeIsNotHiddenBySetupCache(t)
commands = radia.internal.callMex("api.commands");
verifyNotEmpty(t, commands);
verifyError(t, @() radia.setup( ...
    PythonExecutable=fullfile(tempname("C:\temp"), "missing-python.exe"), ...
    ConfigureSimulinkFileGeneration=false), "radia:setup:Python");
verifyEqual(t, radia.internal.callMex("api.commands"), commands);
end

function removeShadow(folder)
if contains(path, folder)
    rmpath(folder);
end
clear radia_mex
if isfolder(folder)
    rmdir(folder, "s");
end
end
