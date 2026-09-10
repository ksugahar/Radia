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
 Duration_s=1e-3,Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,CoreVolume_m3=1e-6,MaxStep_s=5e-6, ...
 MaxIterations=10,RelativeTolerance=1e-2,Relaxation=0.7,CouplingSamples=51,OutputDirectory="C:\temp\radia_hysteretic_interval_test");
verifyTrue(t,r.converged);verifyGreaterThan(t,r.iterations,1);verifyGreaterThan(t,r.current_A(end),0);verifyGreaterThan(t,r.B_T(end),0);verifyGreaterThan(t,max(abs(r.back_emf_V)),0);verifyNotEqual(t,r.hysteresis_state,hys);
verifyEqual(t,r.back_emf_sign_convention,"positive voltage drop in the direction of positive winding current; e = d(N*A*B)/dt");
verifyLessThan(t,abs(r.energy_balance_residual_J)/abs(r.electrical_magnetic_energy_J),0.05);
clear cleanup
end
function testSimulinkPowerElectronicsHysteresisBlock(t)
radia.UtiDelAll();cleanupAll=onCleanup(@()radia.UtiDelAll());fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_hysteretic_drive.cir");
model="radia_power_hysteresis_e2e";cleanupModel=onCleanup(@()closeModel(model));new_system(model);
add_block("simulink/Sources/Constant",model+"/CommandAndPosition",Value="[1;0]",Position=[25 75 105 115]);
radia.simulink.buildHystereticLTspiceBlock(model,Netlist=fixture,Tables={{[0,1,2],[0,1,2]}},EtaOrChi=0.1, ...
 CurrentTrace="I(L1)",Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,CoreVolume_m3=1e-6,SampleTime_s=1e-3, ...
 MaxIterations=10,RelativeTolerance=1e-2,Relaxation=0.7,MaxStep_s=5e-6,CouplingSamples=51,Save=false);
add_block("simulink/Sinks/To Workspace",model+"/CoupledOutputs",VariableName="hys_y",SaveFormat="Array",Position=[410 80 500 115]);
add_line(model,"CommandAndPosition/1","Hysteretic LTspice Plant/1");add_line(model,"Hysteretic LTspice Plant/1","CoupledOutputs/1");set_param(model,StopTime="0.001",Solver="FixedStepDiscrete",FixedStep="0.001");
first=sim(model);second=sim(model);y=first.hys_y;verifySize(t,y,[2,6]);verifyEqual(t,second.hys_y,y,"AbsTol",0);verifyGreaterThan(t,y(1),0);verifyGreaterThan(t,y(2),0);verifyGreaterThan(t,abs(y(4)),0);verifyGreaterThan(t,y(5),0);
logs=dir(fullfile("C:\temp\radia_hysteretic_ltspice_block",'*radia_power_hysteresis_e2e*.log'));verifySize(t,logs,[1,1]);
record=strtrim(string(fileread(fullfile(logs.folder,logs.name))));lines=splitlines(record);
verifyTrue(t,contains(lines(1),"schema=radia.simulink.ltspice_hysteresis.run.v1"));
verifyEqual(t,numel(lines(startsWith(lines,"step="))),2,"One record line per executed step.");
verifyTrue(t,contains(lines(end),"finished failed=0"));
clear cleanupModel cleanupAll
end
function testAirGapMmfAndCouplingSampleConvergence(t)
radia.UtiDelAll();cleanup=onCleanup(@()radia.UtiDelAll());material=radia.MatPlayHysteresis(1,0.1,{{[0,1,2],[0,1,2]}});hys=radia.MatHysSaveState(material);circuit=emptyCircuitState();fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_hysteretic_drive.cir");
samples=[51,101,201];runs=cell(size(samples));
for k=1:numel(samples)
 runs{k}=runReferenceInterval(fixture,material,hys,circuit,samples(k),1e-3);
end
verifyLessThan(t,runs{2}.B_T(end),0.01,"A 1 mm air gap must enter as air reluctance, not as extra iron length.");
mmf=runs{2}.H_A_per_m(end)*1+runs{2}.B_T(end)*1e-3/(4*pi*1e-7);verifyEqual(t,mmf,runs{2}.current_A(end),"RelTol",1e-10);
verifyLessThan(t,abs(runs{2}.B_T(end)-runs{3}.B_T(end)),abs(runs{1}.B_T(end)-runs{2}.B_T(end)));
balance=arrayfun(@(k)abs(runs{k}.energy_balance_residual_J)/abs(runs{k}.electrical_magnetic_energy_J),1:numel(runs));verifyLessThan(t,balance(2),0.01);verifyLessThan(t,balance(3),balance(2));clear cleanup
end
function testEnergyHysteresisIsRejectedAtBlockBoundary(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_hysteretic_drive.cir");model="radia_energy_hysteresis_rejection";cleanup=onCleanup(@()closeModel(model));new_system(model);
verifyError(t,@()radia.simulink.buildHystereticLTspiceBlock(model,Netlist=fixture,Tables={{[0,1],[0,1]}},HysteresisKind="energy",CurrentTrace="I(L1)",Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,CoreVolume_m3=1e-4,SampleTime_s=1e-3,Save=false),"MATLAB:validators:mustBeMember");clear cleanup
end
function state=emptyCircuitState(),state=struct("schema","radia.ltspice.transient_state.v1","time_s",0,"node_names",strings(0,1),"node_voltages_V",zeros(0,1),"inductor_names",strings(0,1),"inductor_currents_A",zeros(0,1));end
function r=runReferenceInterval(fixture,material,hys,circuit,samples,airGap)
r=radia.simulink.runHystereticLTspiceInterval(fixture,material,hys,circuit,CommandValue=1,CurrentTrace="I(L1)",Duration_s=1e-3,Turns=1,CoreArea_m2=1e-6,MagneticPath_m=1,AirGap_m=airGap,CoreVolume_m3=1e-6,MaxStep_s=5e-6,MaxIterations=10,RelativeTolerance=1e-2,Relaxation=0.7,CouplingSamples=samples,OutputDirectory=string(tempname("C:\temp")));
end
function closeModel(name),if bdIsLoaded(name),close_system(name,0);end,end
