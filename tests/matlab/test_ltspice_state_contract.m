function tests = test_ltspice_state_contract
% State extraction contracts that do not require the LTspice executable.
tests = functiontests(localfunctions);
end
function setupOnce(t)
t.TestData.OldPath=path;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),"matlab"));
end
function teardownOnce(t)
path(t.TestData.OldPath);
end
function setup(t)
t.TestData.Folder=string(tempname("C:/temp"));mkdir(t.TestData.Folder);
end
function teardown(t)
rmdir(t.TestData.Folder,'s');
end
function testMissingInternalState(t)
p=netlist(t,["State contract","X1 in 0 dynamic",".subckt dynamic a b", ...
    "L1 a m 1m","R1 m b 1",".ends",".end"]);
verifyError(t,@()radia.ltspice.extractTransientState(raw(),NetlistFile=p), ...
    "radia:ltspice:MissingSubcircuitState");
end
function testStatelessSubcircuit(t)
p=netlist(t,["State contract","X1 in 0 passive params: R=1k", ...
    ".subckt passive a b params: R=1k","R1 a b {R}",".ends",".end"]);
s=radia.ltspice.extractTransientState(raw(),NetlistFile=p);
verifyEmpty(t,s.inductor_names);
end
function testUnsupportedChildDominatesParentCapacitor(t)
p=netlist(t,["State contract","X1 in 0 outer",".subckt outer a b", ...
    "C1 a b 1u","X2 a b inner",".ends",".subckt inner a b", ...
    "D1 a b DM",".model DM D(Cjo=1p)",".ends",".end"]);
verifyError(t,@()radia.ltspice.extractTransientState(raw(),NetlistFile=p), ...
    "radia:ltspice:UnsupportedSubcircuitState");
end
function testHierarchicalInductorNamesSurviveInjection(t)
d=struct("names",["time","V(in)","I(X1:L1)"], ...
    "values",[0,0,0;1e-3,1,2],"step_ranges",[1,2]);
s=radia.ltspice.extractTransientState(d);
p=netlist(t,["State contract","L1 in 0 1m",".end"]);
destination=fullfile(t.TestData.Folder,"state.cir");
radia.ltspice.applyTransientState(p,s,destination,Duration_s=1e-3);
verifyTrue(t,contains(string(fileread(destination)),".ic I(X1:L1)=2"));
end
function d=raw()
d=struct("names",["time","V(in)"],"values",[0,0;1e-3,1],"step_ranges",[1,2]);
end
function p=netlist(t,lines)
p=fullfile(t.TestData.Folder,"source.cir");
fid=fopen(p,'w');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',lines);
end
