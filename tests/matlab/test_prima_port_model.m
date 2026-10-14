function tests = test_prima_port_model
%TEST_PRIMA_PORT_MODEL PRIMA port models from Python into dss and the LTI System block.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeTrue(license("test", "Control_Toolbox") && ~isempty(ver("control")), ...
    "Control System Toolbox is required.");
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
folder = string(tempname("C:\temp"));
pythonExecutable = string(getenv("RADIA_PYTHON_EXECUTABLE"));
if strlength(pythonExecutable) == 0
    pythonExecutable = "python";
end
helper = fullfile(repoRoot, "tests", "matlab", "prima_port_model_python_reference.py");
command = quoteCommandArgument(pythonExecutable) + " " + ...
    quoteCommandArgument(helper) + " " + quoteCommandArgument(folder);
[status, output] = radia.internal.runPythonProcess(command);
if status ~= 0
    error("radia:test:PythonReference", ...
        "Python PRIMA reference failed (%d): %s", status, output);
end
testCase.TestData.folder = folder;
end

function teardownOnce(testCase)
if isfield(testCase.TestData, "folder") && isfolder(testCase.TestData.folder)
    rmdir(testCase.TestData.folder, "s");
end
end

function testLoaderReproducesDescriptorResponse(testCase)
model = radia.simulink.loadPrimaPortModel(fullfile(testCase.TestData.folder, "rl_admittance.json"));
testCase.verifyEqual(model.orientation, "admittance");
testCase.verifyLessThan(model.check_relative_error, 1e-9);
omega = 2 * pi * 700;
direct = model.C * ((1i * omega * model.E - model.A) \ model.B) + model.D;
testCase.verifyEqual(freqresp(model.sys, omega), direct, "RelTol", 1e-12);
testCase.verifyEqual(string(model.sys.InputName), "coil.voltage");
testCase.verifyEqual(string(model.sys.OutputName), "coil.current");
testCase.verifyEqual(string(model.sys.OutputUnit), "A");
end

function testMultiPortImpedanceKeepsPortOrder(testCase)
model = radia.simulink.loadPrimaPortModel(fullfile(testCase.TestData.folder, "rc_impedance.json"));
testCase.verifyEqual(model.port_names, ["node1", "node2"]);
testCase.verifyEqual(string(model.sys.InputName).', ["node1.current", "node2.current"]);
% DC impedance of the RC network is the inverse conductance matrix.
G = -model.A;
testCase.verifyEqual(dcgain(model.sys), inv(G), "RelTol", 1e-12);
end

function testLTIBlockStepResponseMatchesStep(testCase)
model = radia.simulink.loadPrimaPortModel(fullfile(testCase.TestData.folder, "rl_admittance.json"));
name = "prima_lti_" + extractAfter(string(matlab.lang.internal.uuid()), 24);
new_system(name);
testCase.addTeardown(@() close_system(name, 0));
add_block("simulink/Sources/Step", name + "/Step", "Time", "0", "Before", "0", "After", "1");
radia.simulink.addPrimaLTIBlock(name + "/Coil", model);
add_block("simulink/Sinks/Out1", name + "/Current");
add_line(name, "Step/1", "Coil/1");
add_line(name, "Coil/1", "Current/1");
times = [1e-3; 3e-3; 1e-2];
set_param(name, "StopTime", "0.012", "RelTol", "1e-10", "AbsTol", "1e-12", ...
    "SaveOutput", "on", "OutputSaveName", "yout", "SaveFormat", "Dataset", ...
    "OutputOption", "SpecifiedOutputTimes", "OutputTimes", mat2str(times.'));
out = sim(name);
signal = out.yout{1}.Values;
% Exact unit-step response of E x' = A x + B (E invertible here).
F = model.E \ model.A;
G = model.E \ model.B;
expected = arrayfun(@(t) model.C * (F \ ((expm(F * t) - eye(size(F))) * G)), times);
testCase.verifyEqual(signal.Data(ismember(signal.Time, times)), expected, "RelTol", 1e-6);
end

function testImproperOrientationIsRefusedByBlock(testCase)
model = radia.simulink.loadPrimaPortModel(fullfile(testCase.TestData.folder, "rl_impedance.json"));
testCase.verifyFalse(isproper(model.sys));
name = "prima_lti_" + extractAfter(string(matlab.lang.internal.uuid()), 24);
new_system(name);
testCase.addTeardown(@() close_system(name, 0));
testCase.verifyError(@() radia.simulink.addPrimaLTIBlock(name + "/Coil", model), ...
    "radia:simulink:PrimaImproper");
end

function testTamperedCheckResponseIsRejected(testCase)
source = fullfile(testCase.TestData.folder, "rl_admittance.json");
payload = jsondecode(fileread(source));
payload.arrays.A.values(1) = 2 * payload.arrays.A.values(1);
tampered = fullfile(testCase.TestData.folder, "tampered.json");
writelines(jsonencode(payload), tampered);
testCase.verifyError(@() radia.simulink.loadPrimaPortModel(tampered), ...
    "radia:simulink:PrimaExchange");
end

function testInvalidCheckFrequenciesAreRejected(testCase)
for frequency = {[], [100; -1; 1e3], [0; 1e3]}
    payload = localPayload(testCase);
    payload.check.frequency_hz = frequency{1};
    localVerifyRejected(testCase, payload);
end
end

function testInvalidCheckLimitIsRejected(testCase)
for limit = {1e-3, 0, -1e-9, [1e-9, 1e-9]}
    payload = localPayload(testCase);
    payload.check.relative_limit = limit{1};
    localVerifyRejected(testCase, payload);
end
end

function testNonFiniteModelResponseIsRejected(testCase)
% E = A = 0 makes s E - A singular at every check frequency.
payload = localPayload(testCase);
payload.arrays.E.values(:) = 0;
payload.arrays.A.values(:) = 0;
localVerifyRejected(testCase, payload);
end

function payload = localPayload(testCase)
payload = jsondecode(fileread(fullfile(testCase.TestData.folder, "rl_admittance.json")));
end

function localVerifyRejected(testCase, payload)
tampered = string(tempname(testCase.TestData.folder)) + ".json";
writelines(jsonencode(payload), tampered);
testCase.verifyError(@() radia.simulink.loadPrimaPortModel(tampered), ...
    "radia:simulink:PrimaExchange");
end

function value = quoteCommandArgument(value)
value = '"' + replace(string(value), '"', '""') + '"';
end
