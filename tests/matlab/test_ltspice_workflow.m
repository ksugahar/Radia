function tests = test_ltspice_workflow
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(repositoryRoot, "matlab");
entries = string(strsplit(path, pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries, string(matlabDirectory)));
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testInstalledLTspiceRunsAndRawIsParsed(testCase)
fixture = fullfile(fileparts(mfilename("fullpath")), ...
    "fixtures", "ltspice_rc.cir");
result = radia.ltspice.run(fixture, Parameters=struct("Rval", 2000));
verifyEqual(testCase, result.schema, "radia.ltspice.run.v1");
verifyEqual(testCase, result.waveform.names, ...
    ["time", "V(in)", "V(out)", "I(V1)"]);
verifyGreaterThan(testCase, height(result.waveform.values), 100);
verifyTrue(testCase, isfile(result.raw_file));
verifyTrue(testCase, contains(result.log, "Total elapsed time"));
verifyEqual(testCase, result.parameters.Rval, 2000);
end

function testAscSchematicEditConvertAndRun(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.asc");
edited=fullfile("C:\temp","radia_ltspice_edited.asc"); editor=radia.ltspice.SchematicEditor(fixture); editor.setComponentValue("R1","2k"); editor.saveAs(edited);
converted=radia.ltspice.schematicToNetlist(edited,OutputDirectory="C:\temp\radia_ltspice_asc_convert");
verifyTrue(testCase,contains(string(fileread(converted.netlist)),"R1 N001 NC_01 2k"));
result=radia.ltspice.run(edited); verifyEqual(testCase,result.schema,"radia.ltspice.run.v1"); verifyFalse(testCase,isempty(result.schematic_conversion));
end

function testPythonNetlistToSchematicWrapper(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
output=fullfile("C:\temp","radia_ltspice_from_cir.asc");
result=radia.ltspice.netlistToSchematic(fixture,OutputFile=output,ValidateRoundTrip=true);
verifyEqual(testCase,result.schema,"radia.ltspice.netlist_to_schematic.v1");
verifyTrue(testCase,isfile(output)); verifyTrue(testCase,result.validation.topology.equivalent);
verifyTrue(testCase,contains(string(fileread(output)),"SYMBOL"));
end

function testUnknownParameterIsRejected(testCase)
fixture = fullfile(fileparts(mfilename("fullpath")), ...
    "fixtures", "ltspice_rc.cir");
verifyError(testCase, ...
    @() radia.ltspice.run(fixture, Parameters=struct("missing", 1)), ...
    "radia:ltspice:ParameterNotFound");
end

function testPwlExportForSimulinkSignal(testCase)
destination = fullfile("C:\temp", "radia_ltspice_test_gate.pwl");
info = radia.ltspice.writePwl(destination, [0; 1e-6; 2e-6], [0; 1; 0]);
verifyEqual(testCase, info.schema, "radia.ltspice.pwl.v1");
verifyEqual(testCase, info.sample_count, 3);
verifyTrue(testCase, isfile(destination));
text = string(fileread(destination));
verifyTrue(testCase, contains(text, "9.9999999999999995e-07s"));
end

function testSimulinkWaveformDrivesLTspice(testCase)
fixture = fullfile(fileparts(mfilename("fullpath")), ...
    "fixtures", "ltspice_pwl_rc.cir");
gate = [0, 0; 1e-6, 0; 1.01e-6, 1; 10e-6, 1; 10.01e-6, 0; 20e-6, 0];
result = radia.simulink.runLTspice(fixture, ...
    InputSignals=struct("gate", gate));
verifyEqual(testCase, result.schema, "radia.simulink.ltspice.run.v1");
inputIndex = find(result.waveform.names == "V(in)", 1);
verifyGreaterThan(testCase, max(result.waveform.values(:, inputIndex)), 0.99);
verifyEqual(testCase, result.input_pwl.gate.sample_count, height(gate));
end

function testOptunaRunnerEvaluatesCircuitTrial(testCase)
fixture = fullfile(fileparts(mfilename("fullpath")), ...
    "fixtures", "ltspice_rc.cir");
study = radia.optuna.createStudy(direction="maximize", AutoSave=false);
runner = radia.optuna.LTspiceRunner(fixture, ...
    ConfigureFcn=@configureTrial, ScoreFcn=@scoreTrial);
table = runner.optimize(study, 2);
verifyEqual(testCase, table.State, ["COMPLETE"; "COMPLETE"]);
verifyTrue(testCase, all(isfinite(table.Value)));
verifyEqual(testCase, height(study.ParamTable), 2);
end

function testMatlabNativePyLTSpiceEquivalentClasses(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
editor=radia.ltspice.SpiceEditor(fixture); editor.setParameter("Rval",1500);
edited=fullfile("C:\temp","radia_ltspice_edited.cir"); editor.saveAs(edited);
verifyTrue(testCase,contains(string(fileread(edited)),".param Rval=1500"));
runner=radia.ltspice.SimRunner(OutputFolder="C:\temp\radia_ltspice_runner_test");
result=runner.runNow(edited,RunName="single");
verifyTrue(testCase,any(result.raw.getTraceNames()=="V(out)"));
verifyGreaterThan(testCase,numel(result.raw.getTrace("V(out)")),100);
verifyGreaterThan(testCase,result.log_reader.getMeasure("vmax"),0.9);
end

function testBinaryRawAndStepSeparation(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_step_rc.cir");
result=radia.ltspice.run(fixture,RawFormat="binary"); raw=radia.ltspice.RawRead(result.raw_file);
verifyEqual(testCase,raw.Data.schema,"radia.ltspice.raw.binary.v1");
verifyEqual(testCase,raw.getStepCount(),2);
verifyGreaterThan(testCase,numel(raw.getStep("V(out)",1)),100);
verifyGreaterThan(testCase,numel(raw.getStep("V(out)",2)),100);
end

function testSimRunnerMultipleCases(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
runner=radia.ltspice.SimRunner(OutputFolder="C:\temp\radia_ltspice_many_test");
results=runner.runMany(fixture,{struct("Rval",1000),struct("Rval",2000)});
verifyEqual(testCase,numel(results),2); verifyEqual(testCase,results{2}.parameters.Rval,2000);
end

function testComplexAcRawBecomesMatlabComplex(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_ac_rc.cir");
result=radia.ltspice.run(fixture,RawFormat="binary"); raw=radia.ltspice.RawRead(result.raw_file);
frequency=raw.getTrace("frequency"); vout=raw.getTrace("V(out)");
verifyTrue(testCase,raw.isComplex()); verifyEqual(testCase,imag(frequency),zeros(size(frequency)));
expected=1/(1+1i*2*pi*frequency(1)*1000*1e-6);
verifyEqual(testCase,vout(1),expected,"RelTol",1e-12);
verifyLessThan(testCase,angle(vout(end)),0);
verifyEqual(testCase,raw.Data.contract_schema,"radia.ltspice.raw.v2");
verifyEqual(testCase,raw.Data.analysis,"ac");
end

function testComplexAsciiRawIsParsed(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_ac_rc.cir");
result=radia.ltspice.run(fixture,RawFormat="ascii"); raw=radia.ltspice.RawRead(result.raw_file);
verifyTrue(testCase,raw.isComplex()); verifyEqual(testCase,raw.Data.analysis,"ac");
vout=raw.getTrace("V(out)"); verifyLessThan(testCase,angle(vout(end)),0);
end

function testNoiseAndFFTAnalysisAPIs(testCase)
folder=fullfile(fileparts(mfilename("fullpath")),"fixtures");
noise=radia.ltspice.runNoise(fullfile(folder,"ltspice_noise_rc.cir"));
verifyEqual(testCase,noise.schema,"radia.ltspice.noise.v1"); verifyGreaterThan(testCase,numel(noise.frequency_hz),10); verifyFalse(testCase,isempty(fieldnames(noise.noise_traces)));
transient=radia.ltspice.run(fullfile(folder,"ltspice_rc.cir")); fftResult=radia.ltspice.analyzeFFT(transient,"V(out)",SampleCount=1024);
verifyEqual(testCase,fftResult.schema,"radia.ltspice.fft.v1"); verifyEqual(testCase,numel(fftResult.frequency_hz),513); verifyGreaterThanOrEqual(testCase,min(fftResult.amplitude),0);
end

function testRecursiveDependenciesAndStateHandoff(testCase)
fixtureFolder=fullfile(fileparts(mfilename("fullpath")),"fixtures");
fixture=fullfile(fixtureFolder,"ltspice_dependency_root.cir"); manifest=radia.ltspice.collectDependencies(fixture);
verifyEqual(testCase,numel(manifest.local_files),3);
result=radia.ltspice.run(fixture,OutputDirectory="C:\temp\radia_ltspice_dependency_test");
verifyTrue(testCase,isfile(fullfile(result.output_directory,"models","stage1.inc")));
rc=fullfile(fixtureFolder,"ltspice_rc.cir"); intervals=radia.ltspice.runIntervals(rc,[5e-4;5e-4],OutputDirectory="C:\temp\radia_ltspice_interval_test",MaxStep_s=1e-5);
verifyEqual(testCase,intervals.schema,"radia.ltspice.interval_run.v1"); verifyEqual(testCase,numel(intervals.runs),2);
verifyEqual(testCase,intervals.runs{2}.waveform.values(1,1),5e-4,"AbsTol",1e-15);
verifyTrue(testCase,all(ismember(["in";"out"],intervals.states{1}.node_names)));
end

function testParallelOptunaRunnerEvaluatesCircuitTrials(testCase)
if isempty(ver("parallel")), testCase.assumeFail("Parallel Computing Toolbox is unavailable."); end
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
study=radia.optuna.createStudy(AutoSave=false);
runner=radia.optuna.LTspiceRunner(fixture,ConfigureFcn=@configureTrial,ScoreFcn=@scoreTrial);
result=runner.optimizeParallel(study,2,ShowProgress=false);
verifyEqual(testCase,result.State,["COMPLETE";"COMPLETE"]);
end

function parameters = configureTrial(trial)
parameters = struct("Rval", trial.suggestFloat("Rval", 500, 2500));
end
function score = scoreTrial(result, ~)
index = find(result.waveform.names == "V(out)", 1);
score = result.waveform.values(end, index);
end
