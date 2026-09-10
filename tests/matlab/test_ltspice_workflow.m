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
testCase.assumeTrue(radia.ltspice.LTspice.isAvailable(), ...
    "Current ADI LTspice is not installed on this test host.");
testCase.TestData.TempDirectory=string(tempname("C:\temp"));mkdir(testCase.TestData.TempDirectory);
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.MatlabDirectory);
end
if isfolder(testCase.TestData.TempDirectory),rmdir(testCase.TestData.TempDirectory,'s');end
end

function testBinaryRawDoubleFlagIsDecodedWithoutSilentStrideError(testCase)
source=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
netlist=string(fileread(source));netlist=replace(netlist,".tran 0 20m 0 10u",".options numdgt=7"+newline+".tran 0 20m 0 10u");
fixture=fullfile(testCase.TestData.TempDirectory,"double_precision.cir");writeTextFixture(fixture,netlist);
result=radia.ltspice.run(fixture,RawFormat="binary",OutputDirectory=fullfile(testCase.TestData.TempDirectory,"double_run"));
verifyTrue(testCase,any(result.waveform.flags=="double"));
verifyEqual(testCase,result.waveform.values(1,:),zeros(1,4),'AbsTol',0);
verifyGreaterThan(testCase,result.waveform.values(2,1),0);
verifyGreaterThan(testCase,result.waveform.values(2,2),0);
verifyLessThan(testCase,result.waveform.values(2,1),1e-6);
end

function testFastAccessDoubleRawUsesColumnMajorStorage(testCase)
fixture="C:\temp\ltrev\fb.raw";
if ~isfile(fixture),testCase.assumeFail("Opus fast-access fixture is unavailable.");end
raw=radia.ltspice.readRawBinary(fixture);
verifyTrue(testCase,any(raw.flags=="fastaccess"));
verifyEqual(testCase,max(raw.values(:,1)),0.02,'AbsTol',1e-12);
verifyEqual(testCase,max(abs(raw.values(:,2))),1,'AbsTol',1e-12);
end

function testRawPropertiesStopBeforeVariableTable(testCase)
fixture="C:\temp\ltrev\b.raw";
if ~isfile(fixture),testCase.assumeFail("Opus binary fixture is unavailable.");end
raw=radia.ltspice.readRawBinary(fixture);
verifyFalse(testCase,isfield(raw.raw_properties,"Variables"));
verifyFalse(testCase,isfield(raw.raw_properties,"x2I_x1"));
end

function testMissingSubcircuitStateFailsLoudly(testCase)
netlist=fullfile(testCase.TestData.TempDirectory,"sub_state.cir");
writeTextFixture(netlist,"X1 in 0 dynamic"+newline+".subckt dynamic a b"+newline+"L1 a b 1m"+newline+".ends"+newline+".end");
raw=struct("names",["time","V(in)"],"values",[0,0;1e-3,1],"step_ranges",[1,2]);
verifyError(testCase,@()radia.ltspice.extractTransientState(raw,NetlistFile=netlist),"radia:ltspice:MissingSubcircuitState");
end

function testProcessTreeTerminationUsesFrameworkCompatibleApi(testCase)
if ~ispc,testCase.assumeFail("Windows-only process lifecycle test.");end
info=System.Diagnostics.ProcessStartInfo();info.FileName='pwsh';
info.Arguments='-NoLogo -NoProfile -NonInteractive -Command "Start-Sleep -Seconds 30"';
info.UseShellExecute=false;info.CreateNoWindow=true;
process=System.Diagnostics.Process();process.StartInfo=info;verifyTrue(testCase,process.Start());
cleanup=onCleanup(@()radia.ltspice.internal.terminateProcessTree(process));
verifyTrue(testCase,radia.ltspice.internal.terminateProcessTree(process));
verifyTrue(testCase,process.HasExited);clear cleanup
end

function testStateInjectionFailsWithoutEndAndAcceptsHierarchicalInductor(testCase)
fixture=fullfile(testCase.TestData.TempDirectory,"missing_end.cir");writeTextFixture(fixture,"V1 in 0 1"+newline+".tran 1m");
state=struct("schema","radia.ltspice.transient_state.v1","time_s",0,"node_names","in","node_voltages_V",1,"inductor_names","X1:L1","inductor_currents_A",2);
verifyError(testCase,@()radia.ltspice.applyTransientState(fixture,state,fullfile(testCase.TestData.TempDirectory,"out.cir"),Duration_s=1e-3),"radia:ltspice:MissingEnd");
raw=struct("names",["time","V(in)","I(X1:L1)"],"values",[0,0,0;1e-3,1,2],"step_ranges",[1,2]);extracted=radia.ltspice.extractTransientState(raw);
verifyEqual(testCase,extracted.inductor_names,"X1:L1");verifyEqual(testCase,extracted.inductor_currents_A,2);
end

function testIntervalShiftKeepsSignalAndMatrixAxisSynchronized(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
result=radia.ltspice.runIntervals(fixture,[1e-4;1e-4],OutputDirectory=fullfile(testCase.TestData.TempDirectory,"intervals"),MaxStep_s=1e-5);
for k=1:numel(result.runs),verifyEqual(testCase,result.runs{k}.waveform.signals.time,result.runs{k}.waveform.values(:,1),'AbsTol',0);end
end

function testDependenciesSeparateRootSiblingAndExternalFiles(testCase)
root=fullfile(testCase.TestData.TempDirectory,"lib");sibling=fullfile(testCase.TestData.TempDirectory,"lib2");mkdir(root);mkdir(sibling);
external=fullfile(sibling,"external.inc");writeTextFixture(external,".param ext=1");local=fullfile(root,"local.inc");writeTextFixture(local,".param local=1");
main=fullfile(root,"main.cir");writeTextFixture(main,".include local.inc"+newline+".include "+external+newline+".end");manifest=radia.ltspice.collectDependencies(main);
verifyEqual(testCase,numel(manifest.local_files),2);verifyEqual(testCase,numel(manifest.absolute_external),1);verifyEqual(testCase,string(manifest.absolute_external),string(java.io.File(external).getCanonicalPath()));
end

function testOddLengthFftDoublesLastPositiveFrequencyBin(testCase)
n=9;t=(0:n-1)'/n;y=cos(2*pi*4*t);raw=struct("names",["time","V(out)"],"values",[t,y]);result=radia.ltspice.analyzeFFT(raw,"V(out)",SampleCount=n,Window="rectangular");
verifyEqual(testCase,result.amplitude(end),1,'AbsTol',1e-12);
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
edited=tempPath(testCase,"edited.asc"); editor=radia.ltspice.SchematicEditor(fixture); editor.setComponentValue("R1","2k"); editor.saveAs(edited);
converted=radia.ltspice.schematicToNetlist(edited,OutputDirectory=tempPath(testCase,"asc_convert"));
verifyTrue(testCase,contains(string(fileread(converted.netlist)),"R1 N001 NC_01 2k"));
result=radia.ltspice.run(edited); verifyEqual(testCase,result.schema,"radia.ltspice.run.v1"); verifyFalse(testCase,isempty(result.schematic_conversion));
end

function testPythonNetlistToSchematicWrapper(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
output=tempPath(testCase,"from_cir.asc");
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
destination = tempPath(testCase,"gate.pwl");
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
edited=tempPath(testCase,"edited.cir"); editor.saveAs(edited);
verifyTrue(testCase,contains(string(fileread(edited)),".param Rval=1500"));
runner=radia.ltspice.SimRunner(OutputFolder=tempPath(testCase,"runner"));
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
runner=radia.ltspice.SimRunner(OutputFolder=tempPath(testCase,"many"));
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

function testPyLTSpiceRawCompatibilityAndRoundTrip(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_ac_rc.cir");
result=radia.ltspice.run(fixture,RawFormat="ascii"); raw=radia.ltspice.RawRead(result.raw_file);
verifyEqual(testCase,raw.get_trace_names(),raw.getTraceNames());
trace=raw.get_trace("V(out)"); verifyClass(testCase,trace,"radia.ltspice.Trace");
verifyEqual(testCase,trace.get_wave(0),raw.getWave("V(out)",0));
writer=radia.ltspice.RawWrite(); writer.PlotName="AC Analysis";
writer.add_traces_from_raw(raw,{"frequency","V(out)"}); output=tempPath(testCase,"roundtrip.raw"); writer.save(output);
copy=radia.ltspice.RawRead(output); verifyEqual(testCase,copy.getTrace("V(out)"),raw.getTrace("V(out)"),"RelTol",1e-14);
csv=tempPath(testCase,"export.csv"); raw.to_csv(csv,{"frequency","V(out)"},0); verifyTrue(testCase,isfile(csv));
end

function testRawIndicesAreZeroBasedAndAliasesAreSafe(testCase)
fixture=tempPath(testCase,"alias_fixture.raw");
writeAliasRawFixture(fixture);raw=radia.ltspice.RawRead(fixture);
verifyEqual(testCase,raw.getTrace(0),[0;1]);
verifyEqual(testCase,raw.getTrace(1),[1;2]);
verifyEqual(testCase,raw.get_trace(1).Name,"V(in)");
verifyTrue(testCase,any(raw.getTraceNames()=="V(diff)"));
verifyEqual(testCase,raw.getTrace("V(diff)"),[4;8]);
verifyError(testCase,@()raw.getTrace(-1),"radia:ltspice:TraceNotFound");
verifyError(testCase,@()raw.getTrace("V(bad)"), ...
    "radia:ltspice:AliasEvaluation");
writer=radia.ltspice.RawWrite();
writer.addTracesFromRaw(raw,{"time","V(in)"});
verifyEqual(testCase,writer.get_trace(0).Name,"time");
verifyEqual(testCase,writer.get_trace(1).Name,"V(in)");
end

function testUnsupportedSimRunnerCompatibilityOptionsFailLoudly(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");
runner=radia.ltspice.SimRunner();runner.add_command_line_switch("-FastAccess");
verifyError(testCase,@()runner.run_now(fixture), ...
    "radia:ltspice:UnsupportedCompatibility");
verifyError(testCase,@()runner.create_raw_file_with( ...
    "combined.raw",{"V(out)"},[]), ...
    "radia:ltspice:UnsupportedCompatibility");
end

function testPyLTSpiceEditorAliases(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir"); editor=radia.ltspice.SpiceEditor(fixture);
verifyTrue(testCase,any(editor.get_components()=="R1")); verifyEqual(testCase,editor.get_parameter("Rval"),"1000");
editor.set_parameters(struct("Rval",2200)); verifyEqual(testCase,editor.get_parameter("Rval"),"2200");
verifyEqual(testCase,editor.get_component_nodes("R1"),["in";"out"]);
end

function testAscGraphicalEditingCompatibility(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.asc");editor=radia.ltspice.AscEditor(fixture);
[position,rotation]=editor.get_component_position("R1");verifyEqual(testCase,position,[160,80]);verifyEqual(testCase,rotation,"R90");
editor.set_component_position("R1",[192,112],"R0");editor.set_component_attribute("R1","SpiceLine","temp=25");editor.addWire([0,0],[16,0]);editor.set_parameter("gain",2);
output=tempPath(testCase,"asc_editor_compat.asc");editor.save_as(output);text=string(fileread(output));
verifyTrue(testCase,contains(text,"SYMBOL res 192 112 R0"));verifyTrue(testCase,contains(text,"SYMATTR SpiceLine temp=25"));verifyTrue(testCase,contains(text,"WIRE 0 0 16 0"));verifyEqual(testCase,editor.get_parameter("gain"),"2");
end

function testSteppedLogQueriesAndRawStepConditions(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_step_rc.cir");result=radia.ltspice.run(fixture,RawFormat="binary");
log=radia.ltspice.LTSpiceLogReader(result.log_file);verifyTrue(testCase,log.has_steps());verifyEqual(testCase,log.get_step_vars(),"rval");verifyEqual(testCase,log.steps_with_parameter_equal_to("rval",2000),1);
raw=radia.ltspice.RawRead(result.raw_file);verifyEqual(testCase,raw.get_steps(struct("rval",1000)),0);
end

function testAsynchronousSimRunnerTask(testCase)
fixture=fullfile(fileparts(mfilename("fullpath")),"fixtures","ltspice_rc.cir");runner=radia.ltspice.SimRunner(OutputFolder=tempPath(testCase,"async"));
task=runner.run(fixture);verifyClass(testCase,task,"radia.ltspice.RunTask");verifyTrue(testCase,task.wait(30));files=task.wait_results();verifyTrue(testCase,isfile(files{1}));verifyTrue(testCase,isfile(files{2}));verifyEqual(testCase,task.Status,"completed");
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
result=radia.ltspice.run(fixture,OutputDirectory=tempPath(testCase,"dependency"));
verifyTrue(testCase,isfile(fullfile(result.output_directory,"models","stage1.inc")));
rc=fullfile(fixtureFolder,"ltspice_rc.cir"); intervals=radia.ltspice.runIntervals(rc,[5e-4;5e-4],OutputDirectory=tempPath(testCase,"interval"),MaxStep_s=1e-5);
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

function writeTextFixture(path,text)
handle=fopen(path,'w');assert(handle>=0);cleanup=onCleanup(@()fclose(handle));fprintf(handle,'%s',text);clear cleanup
end

function path=tempPath(testCase,name),path=fullfile(testCase.TestData.TempDirectory,name);end


function writeAliasRawFixture(path)
lines=[ ...
    "Title: Radia alias regression"; ...
    "Date: 2026-08-04"; ...
    "Plotname: Transient Analysis"; ...
    "Flags: real"; ...
    "Alias: V(diff)=2*(V(out)-V(in))"; ...
    "Alias: V(bad)=V(out);1"; ...
    "No. Variables: 3"; ...
    "No. Points: 2"; ...
    "Variables:"; ...
    "0 time time"; ...
    "1 V(in) voltage"; ...
    "2 V(out) voltage"; ...
    "Values:"; ...
    "0 0"; ...
    "1"; ...
    "3"; ...
    "1 1"; ...
    "2"; ...
    "6"];
writelines(lines,path,Encoding="UTF-8");
end
