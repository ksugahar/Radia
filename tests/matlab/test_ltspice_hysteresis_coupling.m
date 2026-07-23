function tests=test_ltspice_hysteresis_coupling
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
end
function teardownOnce(t)
if t.TestData.RemoveMatlabDirectory
    rmpath(fullfile(t.TestData.Root,"matlab"));
end
end
function testWaveformRelaxedCircuitHysteresisInterval(t)
radia.UtiDelAll();cleanup=onCleanup(@()radia.UtiDelAll());material=radia.MatPlayHysteresis(1,0.1,{{[0,1,2],[0,1,2]}});hys=radia.MatHysSaveState(material);
circuit=emptyCircuitState();fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_hysteretic_drive.cir");
r=radia.simulink.runHystereticLTspiceInterval(fixture,material,hys,circuit,CommandValue=1,CurrentTrace="I(L1)", ...
 Duration_s=1e-3,Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,CoreVolume_m3=1e-4,MaxStep_s=5e-6, ...
 MaxIterations=10,RelativeTolerance=1e-2,Relaxation=0.7,CouplingSamples=51,OutputDirectory="C:\temp\radia_hysteretic_interval_test");
verifyTrue(t,r.converged);verifyGreaterThan(t,r.iterations,1);verifyGreaterThan(t,r.current_A(end),0);verifyGreaterThan(t,r.B_T(end),0);verifyGreaterThan(t,max(abs(r.back_emf_V)),0);verifyNotEqual(t,r.hysteresis_state,hys);
clear cleanup
end
function testSimulinkPowerElectronicsHysteresisBlock(t)
radia.UtiDelAll();cleanupAll=onCleanup(@()radia.UtiDelAll());fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_hysteretic_drive.cir");
model="radia_power_hysteresis_e2e";cleanupModel=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/CommandAndPosition",Value="[1;0]",Position=[25 75 105 115]);
radia.simulink.buildHystereticLTspiceBlock(model,Netlist=fixture,Tables={{[0,1,2],[0,1,2]}},EtaOrChi=0.1, ...
 CurrentTrace="I(L1)",Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,CoreVolume_m3=1e-4,SampleTime_s=1e-3, ...
 MaxIterations=10,RelativeTolerance=1e-2,Relaxation=0.7,MaxStep_s=5e-6,CouplingSamples=51,Save=false);
add_block("simulink/Sinks/To Workspace",model+"/CoupledOutputs",VariableName="hys_y",SaveFormat="Array",Position=[410 80 500 115]);
add_line(model,"CommandAndPosition/1","Hysteretic LTspice Plant/1");add_line(model,"Hysteretic LTspice Plant/1","CoupledOutputs/1");set_param(model,StopTime="0",Solver="FixedStepDiscrete",FixedStep="0.001");
out=sim(model);y=out.hys_y;verifySize(t,y,[1,6]);verifyGreaterThan(t,y(1),0);verifyGreaterThan(t,y(2),0);verifyGreaterThan(t,abs(y(4)),0);verifyGreaterThan(t,y(5),0);
clear cleanupModel cleanupAll
end
function state=emptyCircuitState(),state=struct("schema","radia.ltspice.transient_state.v1","time_s",0,"node_names",strings(0,1),"node_voltages_V",zeros(0,1),"inductor_names",strings(0,1),"inductor_currents_A",zeros(0,1));end
function closeModel(name),if bdIsLoaded(name),close_system(name,0);end,end
