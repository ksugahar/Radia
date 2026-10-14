function tests = test_ltspice_rlc_resonance_optimization
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.netlist = fullfile(repoRoot, "matlab", "samples", ...
    "ltspice_rlc_resonance_ac.cir");
testCase.TestData.outputRoot = "C:\temp\radia_ltspice_rlc_test";
end

function testACResonanceMatchesAnalyticIdentity(testCase)
resistance = 2;
inductance = 10e-3;
target_Hz = 1000;
capacitance = 1/(inductance*(2*pi*target_Hz)^2);
result = radia.ltspice.run(testCase.TestData.netlist, ...
    Parameters=struct("Rval", resistance, "Lval", inductance, ...
        "Cval", capacitance), ...
    OutputDirectory=fullfile(testCase.TestData.outputRoot, "identity"), ...
    RawFormat="binary");
[measured_Hz, details] = radia.ltspice.measureSeriesResonance(result);
verifyEqual(testCase, measured_Hz, target_Hz, "RelTol", 5e-5);
verifyLessThan(testCase, details.zero_vs_peak_relative_difference, 5e-5);
verifyEqual(testCase, details.input_resistance_at_resonance_Ohm, ...
    resistance, "RelTol", 1e-10);
end

function testObjectiveRecordsComplexACMeasurement(testCase)
study = radia.optuna.createStudy( ...
    direction="minimize", sampler=radia.optuna.RandomSampler(7), ...
    AutoSave=false);
trial = study.ask();
value = radia.ltspice.rlcResonanceObjective( ...
    trial, testCase.TestData.netlist, ...
    OutputRoot=fullfile(testCase.TestData.outputRoot, "objective"));
verifyTrue(testCase, isfinite(value));
verifyTrue(testCase, isfield(trial.UserAttrs, "ltspice_ac_resonance"));
measurement = trial.UserAttrs.ltspice_ac_resonance;
verifyGreaterThan(testCase, measurement.frequency_Hz, 600);
verifyLessThan(testCase, measurement.frequency_Hz, 1700);
study.tell(trial, value);
verifyEqual(testCase, study.TrialTable.State, "COMPLETE");
end
