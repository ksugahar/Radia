function tests = test_rlc_resonance_optimization_demo
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.repoRoot = repoRoot;
testCase.TestData.outputDirectory = "C:\temp\radia_rlc_resonance_test";
if isfolder(testCase.TestData.outputDirectory)
    rmdir(testCase.TestData.outputDirectory, "s");
end
mkdir(testCase.TestData.outputDirectory);
end

function teardownOnce(testCase)
closeLoaded("radia_rlc_resonance_test_plant");
if isfolder(testCase.TestData.outputDirectory)
    rmdir(testCase.TestData.outputDirectory, "s");
end
end

function testAnalyticDefinition(testCase)
result = radia.simulink.analyzeRLCResonance(2, 10e-3, 2.5e-6);
verifyGreaterThan(testCase, result.natural_frequency_Hz, ...
    result.ringdown_frequency_Hz);
verifyEqual(testCase, result.series_current_peak_frequency_Hz, ...
    result.natural_frequency_Hz, "AbsTol", eps(result.natural_frequency_Hz));
verifyGreaterThan(testCase, result.quality_factor, 10);
end

function testZeroCrossingMeasurement(testCase)
time_s = (0:2e-6:0.012).';
expected_Hz = 1000;
values = exp(-100*time_s) .* cos(2*pi*expected_Hz*time_s);
[measured_Hz, details] = radia.simulink.measureRingdownFrequency( ...
    [time_s values]);
verifyEqual(testCase, measured_Hz, expected_Hz, "RelTol", 1e-7);
verifyGreaterThan(testCase, details.crossing_count, 20);
end

function testPlantMatchesTheory(testCase)
modelFile = radia.simulink.buildRLCResonanceOptimizationDemo( ...
    OutputDirectory=testCase.TestData.outputDirectory, ...
    ModelName="radia_rlc_resonance_test_plant");
verifyTrue(testCase, isfile(modelFile));

resistance = 2;
inductance = 10e-3;
capacitance = 2.5324e-6;
expected = radia.simulink.analyzeRLCResonance( ...
    resistance, inductance, capacitance);
in = Simulink.SimulationInput(modelFile);
in = in.setVariable("R_ohm", resistance);
in = in.setVariable("L_H", inductance);
in = in.setVariable("C_F", capacitance);
out = sim(in);
[measured_Hz, details] = radia.simulink.measureRingdownFrequency( ...
    out.get("ring_voltage"));
verifyEqual(testCase, measured_Hz, expected.ringdown_frequency_Hz, ...
    "RelTol", 1e-4);
verifyLessThan(testCase, details.relative_half_period_mad, 1e-5);

load_system(modelFile);
verifyEqual(testCase, string(get_param( ...
    "radia_rlc_resonance_test_plant/Series RLC", "BlockType")), ...
    "StateSpace");
verifyEqual(testCase, string(get_param( ...
    "radia_rlc_resonance_test_plant/Ring voltage", "VariableName")), ...
    "ring_voltage");
close_system("radia_rlc_resonance_test_plant", 0);
end

function closeLoaded(model)
if bdIsLoaded(model)
    close_system(model, 0);
end
end
