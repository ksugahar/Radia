function tests=test_ltspice_simulink_coupling
tests=functiontests(localfunctions);
end
function setupOnce(t)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));addpath(fullfile(root,"matlab"));t.TestData.Root=root;
end
function teardownOnce(t),rmpath(fullfile(t.TestData.Root,"matlab"));end

function testStatefulMultipleIOAndRestartCleanup(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");
model="radia_ltspice_stateful_e2e";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Inputs",Value="[1;2]",Position=[30 70 80 110]);
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"], ...
 OutputTraces=["V(out1)";"V(out2)"],SampleTime_s=1e-3,MaxStep_s=1e-5,Timeout_s=30,Save=false);
set_param(model+"/LTspice Circuit","Position",[130 65 300 115]);
add_block("simulink/Sinks/To Workspace",model+"/Results",VariableName="coupled_y",SaveFormat="Array",Position=[350 70 440 110]);
add_line(model,"Inputs/1","LTspice Circuit/1");add_line(model,"LTspice Circuit/1","Results/1");set_param(model,StopTime="0.002",Solver="FixedStepDiscrete",FixedStep="0.001");
first=sim(model);second=sim(model);a=first.coupled_y;b=second.coupled_y;
verifyEqual(t,size(a,2),2);verifyEqual(t,a,b,"AbsTol",1e-12);verifyGreaterThan(t,a(end,1),a(1,1));verifyGreaterThan(t,a(end,2),a(1,2));verifyGreaterThan(t,a(end,2),a(end,1));
clear cleanup
end

function testMissingTraceBecomesSimulinkDiagnostic(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_multi_io.cir");model="radia_ltspice_error_e2e";cleanup=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/Input",Value="[1;2]");
radia.simulink.buildLTspiceBlock(model,Netlist=fixture,InputNames=["drive1";"drive2"],OutputTraces="V(missing)",SampleTime_s=1e-3,Save=false);
add_line(model,"Input/1","LTspice Circuit/1");set_param(model,StopTime="0",Solver="FixedStepDiscrete",FixedStep="0.001");
verifyError(t,@()sim(model),"Simulink:blocks:MSFB_BlockMethodFailed");clear cleanup
end
function closeModel(name),if bdIsLoaded(name),close_system(name,0);end,end
