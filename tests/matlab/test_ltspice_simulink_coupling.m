function tests=test_ltspice_simulink_coupling
tests=functiontests(localfunctions);
end
function setupOnce(t)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
t.TestData.RemoveMatlabDirectory= ...
    ~any(strcmpi(entries,string(matlabDirectory)));
if t.TestData.RemoveMatlabDirectory,addpath(matlabDirectory);end
t.TestData.Root=root;
t.assumeTrue(radia.ltspice.LTspice.isAvailable(),"Current ADI LTspice is not installed on this test host.");
t.TestData.TempDirectory=string(tempname("C:\temp"));mkdir(t.TestData.TempDirectory);
end
function teardownOnce(t)
if t.TestData.RemoveMatlabDirectory
    rmpath(fullfile(t.TestData.Root,"matlab"));
end
if isfolder(t.TestData.TempDirectory),rmdir(t.TestData.TempDirectory,'s');end
end

function testStatefulMultipleIOAndRestartCleanup(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");
model="radia_ltspice_stateful_e2e";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Inputs",Value="[1;2]",Position=[30 70 80 110]);
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"], ...
 OutputTraces=["V(out1)";"V(out2)"],SampleTime_s=1e-3,MaxStep_s=1e-5,Timeout_s=30,RootDirectory=t.TestData.TempDirectory,Save=false);
set_param(model+"/LTspice Circuit","Position",[130 65 300 115]);
add_block("simulink/Sinks/To Workspace",model+"/Results",VariableName="coupled_y",SaveFormat="Array",Position=[350 70 440 110]);
add_line(model,"Inputs/1","LTspice Circuit/1");add_line(model,"LTspice Circuit/1","Results/1");set_param(model,StopTime="0.002",Solver="FixedStepDiscrete",FixedStep="0.001");
first=sim(model);second=sim(model);a=first.coupled_y;b=second.coupled_y;
verifyEqual(t,size(a,2),2);verifyEqual(t,a,b,"AbsTol",1e-12);verifyGreaterThan(t,a(end,1),a(1,1));verifyGreaterThan(t,a(end,2),a(1,2));verifyGreaterThan(t,a(end,2),a(end,1));
verifyEmpty(t,dir(fullfile(t.TestData.TempDirectory,'tp*')));
clear cleanup
end

function testMissingTraceBecomesSimulinkDiagnostic(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");model="radia_ltspice_error_e2e";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Input",Value="[1;2]");
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"],OutputTraces="V(missing)",SampleTime_s=1e-3,RootDirectory=t.TestData.TempDirectory,Save=false);
add_line(model,"Input/1","LTspice Circuit/1");set_param(model,StopTime="0",Solver="FixedStepDiscrete",FixedStep="0.001");
verifyError(t,@()sim(model),"Simulink:blocks:MSFB_BlockMethodFailed");clear cleanup
end

function testFailedRunKeepsArtifactsAndRecordsTheFailure(t)
% A failed run must stay inspectable: the policy requires a log regardless of
% exit code, and the step folder that produced the failure must survive.
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");
root=fullfile(t.TestData.TempDirectory,"failed_run");mkdir(root);
model="radia_ltspice_failed_artifacts";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Input",Value="[1;2]");
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"],OutputTraces="V(missing)",SampleTime_s=1e-3,RootDirectory=root,Save=false);
add_line(model,"Input/1","LTspice Circuit/1");set_param(model,StopTime="0",Solver="FixedStepDiscrete",FixedStep="0.001");
verifyError(t,@()sim(model),"Simulink:blocks:MSFB_BlockMethodFailed");
logs=dir(fullfile(root,'*.log'));verifySize(t,logs,[1,1]);
text=string(fileread(fullfile(logs.folder,logs.name)));
verifyTrue(t,contains(text,"schema=radia.simulink.ltspice.run.v1"));
verifyTrue(t,contains(text,"FAILED"));
verifyTrue(t,contains(text,"finished failed=1"));
verifyNotEmpty(t,dir(fullfile(root,'tp*')));
clear cleanup
end

function testSuccessfulRunKeepsOnlyTheBoundedRecord(t)
% A completed run keeps the small record but not the heavy step folders, and
% the per-step artifact ring stays bounded while the run advances.
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");
root=fullfile(t.TestData.TempDirectory,"bounded_run");mkdir(root);
model="radia_ltspice_bounded_artifacts";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Inputs",Value="[1;2]",Position=[30 70 80 110]);
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"], ...
 OutputTraces=["V(out1)";"V(out2)"],SampleTime_s=1e-3,MaxStep_s=1e-5,Timeout_s=30,RootDirectory=root,Save=false);
set_param(model+"/LTspice Circuit","Position",[130 65 300 115]);
add_line(model,"Inputs/1","LTspice Circuit/1");set_param(model,StopTime="0.005",Solver="FixedStepDiscrete",FixedStep="0.001");
sim(model);
logs=dir(fullfile(root,'*.log'));verifySize(t,logs,[1,1]);
lines=splitlines(strtrim(string(fileread(fullfile(logs.folder,logs.name)))));
steps=lines(startsWith(lines,"step="));
verifyEqual(t,numel(steps),6,"One record line per executed step.");
verifyTrue(t,contains(lines(end),"finished failed=0"));
verifyEmpty(t,dir(fullfile(root,'tp*')),"A completed run must not leave step artifacts behind.");
clear cleanup
end
function testStepArtifactRingKeepsOnlyTheMostRecentFolders(t)
% The ring is what bounds peak disk usage while a long run is still going.
root=fullfile(t.TestData.TempDirectory,"ring");mkdir(root);
kept=strings(0,1);
for k=1:7
    folder=fullfile(root,sprintf("step_%06d",k));mkdir(folder);
    kept(end+1,1)=string(folder); %#ok<AGROW>
    kept=radia.simulink.internal.runArtifacts("prune",kept,3);
    verifyLessThanOrEqual(t,numel(kept),3);
end
verifyEqual(t,numel(kept),3);
verifyTrue(t,all(isfolder(kept)),"The retained folders must still exist.");
remaining=string({dir(fullfile(root,'step_*')).name}');
verifyEqual(t,sort(remaining),sort(["step_000005";"step_000006";"step_000007"]));
end

function closeModel(name),if bdIsLoaded(name),close_system(name,0);end,end
