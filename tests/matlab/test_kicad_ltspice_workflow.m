function tests=test_kicad_ltspice_workflow
tests=functiontests(localfunctions);
end
function setupOnce(t),t.TestData.Root=fileparts(fileparts(fileparts(mfilename("fullpath"))));end
function testSpiceExportContract(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","kicad_ltspice_bridge.kicad_sch");
cli=fullfile(t.TestData.Root,"tests","matlab","fixtures","fake_kicad_cli.cmd");
output=fullfile("C:\temp","radia_kicad_test","bridge.cir");
r=radia.kicad.exportSpiceNetlist(fixture,OutputFile=output,Executable=cli);
verifyEqual(t,r.schema,"radia.kicad.spice_netlist.v1");verifyEqual(t,r.format,"spice");
verifyTrue(t,isfile(r.netlist_file));verifySubstring(t,fileread(r.netlist_file),"R1 in out 1k");
end
function testRejectsNonKiCadSchematic(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","ltspice_rc.cir");
verifyError(t,@()radia.kicad.exportSpiceNetlist(fixture),"radia:kicad:SchematicRequired");
end
function testPrepareLTspiceRoundTrip(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","kicad_ltspice_bridge.kicad_sch");
cli=fullfile(t.TestData.Root,"tests","matlab","fixtures","fake_kicad_cli.cmd");
r=radia.kicad.prepareLTspice(fixture,OutputDirectory=fullfile("C:\temp","radia_kicad_roundtrip"), ...
 KiCadExecutable=cli,ValidateRoundTrip=true);
verifyEqual(t,r.schema,"radia.kicad.ltspice_preparation.v1");
verifyTrue(t,isfile(r.netlist_file));verifyTrue(t,isfile(r.ltspice_schematic));
verifyTrue(t,r.ltspice_conversion.validation.topology.equivalent);
end
function testBuildsSimulinkLTspiceBlock(t)
fixture=fullfile(t.TestData.Root,"tests","matlab","fixtures","kicad_ltspice_bridge.kicad_sch");
cli=fullfile(t.TestData.Root,"tests","matlab","fixtures","fake_kicad_cli.cmd");model="radia_kicad_ltspice_test";
cleanup=onCleanup(@()closeModel(model));
r=radia.kicad.buildLTspiceBlock(model,fixture,OutputDirectory=fullfile("C:\temp","radia_kicad_block"), ...
 KiCadExecutable=cli,InputNames="drive",OutputTraces="V(out)",SampleTime_s=1e-3,Save=false);
verifyEqual(t,r.schema,"radia.kicad.simulink_ltspice.v1");
verifyEqual(t,string(get_param(r.block_path,"FunctionName")),"radia_ltspice_sfun");
verifySubstring(t,string(get_param(r.block_path,"Parameters")),string(r.netlist_file));
clear cleanup
end
function closeModel(name),if bdIsLoaded(name),close_system(name,0);end,end
