function tests=test_ltspice_data_safety
tests=functiontests(localfunctions);
end
function setupOnce(t)
t.TestData.OldPath=path;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),"matlab"));
end
function teardownOnce(t),path(t.TestData.OldPath);end
function setup(t)
t.TestData.Folder=string(tempname("C:/temp")); mkdir(t.TestData.Folder);
end
function teardown(t),rmdir(t.TestData.Folder,'s');end
function testDouble(t)
v=[0 1e-100 -3;1e-10 pi 2;0.009792 -1e30 4];
p=fixture(t,"real double forward",v,false);
d=radia.ltspice.readRaw(p);
verifyEqual(t,d.values,v); verifyEqual(t,d.raw_properties.Flags,"real double forward");
verifyEqual(t,d.raw_properties.Offset,"0");
end
function testRealLTspiceDoubleAgainstAscii(t)
folder=fullfile(fileparts(mfilename('fullpath')),"fixtures");
b=radia.ltspice.readRaw(fullfile(folder,"ltspice_numdgt7_binary.raw"));
a=radia.ltspice.readRaw(fullfile(folder,"ltspice_numdgt7_ascii.raw"));
verifyTrue(t,contains(b.raw_properties.Flags,"double"));
verifyEqual(t,b.names,a.names);
verifyEqual(t,b.values,a.values,RelTol=1e-12,AbsTol=1e-15);
end
function testMixedPrecision(t)
v=[0 1.25 -3;1e-10 2.5 4];
d=radia.ltspice.readRaw(fixture(t,"real forward",v,false));
verifyEqual(t,d.values,v);
end
function testComplex(t)
v=[1+0i 2+3i;4+0i -5-6i];
d=radia.ltspice.readRaw(fixture(t,"complex forward",v,false));
verifyEqual(t,d.values,v);
end
function testCrLf(t)
v=[0 1;1 2];
d=radia.ltspice.readRaw(fixture(t,"real double",v,true));
verifyEqual(t,d.values,v);
end
function testUnsupportedFlags(t)
p=fixture(t,"real fastaccess",[0 1;1 2],false);
verifyError(t,@()radia.ltspice.readRaw(p),"radia:ltspice:RawFlags");
end
function testConflictingFlags(t)
p=fixture(t,"real complex",[0 1;1 2],false);
verifyError(t,@()radia.ltspice.readRaw(p),"radia:ltspice:RawFlags");
end
function testExtraPayload(t)
p=fixture(t,"real double",[0 1;1 2],false);
f=fopen(p,'ab');fwrite(f,uint8(1),'uint8');fclose(f);
verifyError(t,@()radia.ltspice.readRaw(p),"radia:ltspice:RawFormat");
end
function testTruncated(t)
p=fixture(t,"real double",[0 1;1 2],false);
f=fopen(p,'rb');b=fread(f,inf,'*uint8');fclose(f);
f=fopen(p,'wb');fwrite(f,b(1:end-1),'uint8');fclose(f);
verifyError(t,@()radia.ltspice.readRaw(p),"radia:ltspice:RawTruncated");
end
function testHierarchicalStateFails(t)
d=struct('names',["time","V(a)","I(X1:L1)"], ...
 'values',[0 1 2;1 3 4],'step_ranges',[1 2]);
verifyError(t,@()radia.ltspice.extractTransientState(d), ...
 "radia:ltspice:UnsupportedHierarchicalState");
end
function testTopLevelState(t)
d=struct('names',["time","V(a)","I(L1)","I(V1)"], ...
 'values',[0 1 2 9;1 3 4 9],'step_ranges',[1 2]);
s=radia.ltspice.extractTransientState(d);
verifyEqual(t,s.inductor_names,"L1");verifyEqual(t,s.inductor_currents_A,4);
end
function testMissingEnd(t)
p=fullfile(t.TestData.Folder,"a.cir");writeText(p,"* test"+newline+"R1 a 0 1");
verifyError(t,@()radia.ltspice.applyTransientState(p,state(), ...
 fullfile(t.TestData.Folder,"out.cir"),Duration_s=1),"radia:ltspice:StateNetlistEnd");
end
function testInjection(t)
p=fullfile(t.TestData.Folder,"a.cir");
writeText(p,"* test"+newline+".ic V(a)=999"+newline+"  .end  "+newline);
out=fullfile(t.TestData.Folder,"out.cir");
radia.ltspice.applyTransientState(p,state(),out,Duration_s=1);
text=string(fileread(out));verifyTrue(t,contains(text,".ic V(a)=2"));
verifyFalse(t,contains(text,"999"));verifyTrue(t,contains(text,".tran 0 1 uic"));
end
function s=state()
s=struct('schema',"radia.ltspice.transient_state.v1",'time_s',1, ...
 'node_names',"a",'node_voltages_V',2,'inductor_names',strings(0,1), ...
 'inductor_currents_A',zeros(0,1));
end
function writeText(p,text)
f=fopen(p,'w');c=onCleanup(@()fclose(f));fprintf(f,'%s',text);
end
function p=fixture(t,flags,v,crlf)
p=fullfile(t.TestData.Folder,"test.raw");
eol=string(newline);if crlf,eol=string(char(13))+newline;end
header="Title: deterministic layout regression"+eol+"Flags: "+flags+eol+ ...
 "No. Variables: "+size(v,2)+eol+"No. Points: "+size(v,1)+eol+ ...
 "Offset: 0"+eol+"Variables:"+eol+"0 time time"+eol;
for k=2:size(v,2),header=header+(k-1)+" V(n"+k+") voltage"+eol;end
header=header+"Binary:"+eol;
f=fopen(p,'wb','ieee-le');c=onCleanup(@()fclose(f));
fwrite(f,unicode2native(header,'UTF-16LE'),'uint8');
for row=1:size(v,1)
 if contains(flags,"complex")
  for k=1:size(v,2),fwrite(f,[real(v(row,k)) imag(v(row,k))],'double');end
 elseif contains(flags,"double"),fwrite(f,v(row,:),'double');
 else,fwrite(f,v(row,1),'double');fwrite(f,v(row,2:end),'single');end
end
end
