function tests = test_winding_dictionary_simulink
%TEST_WINDING_DICTIONARY_SIMULINK Circuit and mechanical coupling contracts.
tests = functiontests(localfunctions);
end

function testCompilesSignedCoilSidesAndTerminals(testCase)
materials = materialContract();
phaseA = radia.simulink.makeWindingSpec( ...
    Regions=["phase_a_pos";"phase_a_neg"],RegionPolarity=[1;-1], ...
    Turns=120,ParallelPaths=2,Resistance_ohm=0.4, ...
    PositiveTerminal="a",NegativeTerminal="neutral");
phaseB = radia.simulink.makeWindingSpec( ...
    Regions=["phase_b_pos";"phase_b_neg"],RegionPolarity=[1;-1], ...
    Turns=120,Polarity=-1,Resistance_ohm=0.4, ...
    PositiveTerminal="b",NegativeTerminal="neutral");
windings = dictionary(["phase_a","phase_b"],{phaseA,phaseB});
contract = radia.simulink.compileWindingDictionary(windings,materials, ...
    MaxWindings=4,MaxRegionsPerWinding=4,MaxTerminals=4);
verifyEqual(testCase,contract.schema,"radia.simulink.winding-dictionary.v1");
verifyEqual(testCase,contract.winding_names,["phase_a";"phase_b"]);
verifyEqual(testCase,contract.terminal_names,["a";"b";"neutral"]);
verifyEqual(testCase,contract.runtime.region_id(1:2,1),uint32([3;2]));
verifyEqual(testCase,contract.runtime.region_polarity(1:2,1),int8([1;-1]));
verifyEqual(testCase,contract.runtime.effective_turns(1),60);
verifyEqual(testCase,contract.runtime.effective_turns(2),-120);
verifyFalse(testCase,contract.runtime_policy.dictionary_lookup_per_step);
verifyEqual(testCase,contract.sign_convention.positive_torque,"increasing rotor angle");
bus = radia.simulink.makeWindingBusObject(contract.runtime,AssignToBase=false);
verifyEqual(testCase,numel(bus.Elements),numel(fieldnames(contract.runtime)));
end

function testRejectsUnknownAndSharedRegions(testCase)
materials = materialContract();
unknown = dictionary("phase",{radia.simulink.makeWindingSpec( ...
    Regions="missing",Turns=10)});
verifyError(testCase,@() radia.simulink.compileWindingDictionary(unknown,materials), ...
    "radia:simulink:WindingUnknownRegion");
one = radia.simulink.makeWindingSpec(Regions="phase_a_pos",Turns=10, ...
    PositiveTerminal="a",NegativeTerminal="n");
two = radia.simulink.makeWindingSpec(Regions="phase_a_pos",Turns=10, ...
    PositiveTerminal="b",NegativeTerminal="n");
shared = dictionary(["a","b"],{one,two});
verifyError(testCase,@() radia.simulink.compileWindingDictionary(shared,materials), ...
    "radia:simulink:WindingSharedRegion");
end

function testRejectsBadPolarityAndCapacity(testCase)
verifyError(testCase,@() radia.simulink.makeWindingSpec( ...
    Regions=["one";"two"],RegionPolarity=[1;0],Turns=10), ...
    "radia:simulink:WindingRegionPolarity");
materials = materialContract();
winding = dictionary("phase",{radia.simulink.makeWindingSpec( ...
    Regions=["phase_a_pos";"phase_a_neg"],Turns=10)});
verifyError(testCase,@() radia.simulink.compileWindingDictionary( ...
    winding,materials,MaxRegionsPerWinding=1), ...
    "radia:simulink:WindingRegionCapacity");
end

function testElectromechanicalBusesAreFixedAndSolverIndependent(testCase)
buses = radia.simulink.makeElectromechanicalBusObjects( ...
    MaxWindings=6,AssignToBase=false);
verifyEqual(testCase,size(buses.command_value.terminal_voltage_V),[6 1]);
verifyEqual(testCase,size(buses.response_value.flux_linkage_Wb_turn),[6 1]);
verifyEqual(testCase,size(buses.command_value.translation_position_m),[3 1]);
verifyEqual(testCase,size(buses.response_value.electromagnetic_force_N),[3 1]);
verifyTrue(testCase,isa(buses.command,"Simulink.Bus"));
verifyTrue(testCase,isa(buses.response,"Simulink.Bus"));
end

function testRegistersWindingAndDynamicBuses(testCase)
cleanup = onCleanup(@cleanupWorkspace); %#ok<NASGU>
materials = materialContract();
winding = dictionary("phase",{radia.simulink.makeWindingSpec( ...
    Regions=["phase_a_pos";"phase_a_neg"],RegionPolarity=[1;-1],Turns=10)});
assignin("base","radia_windings_qa",winding);
assignin("base","radia_material_contract_qa",materials);
contract = radia.simulink.registerWindingDictionary( ...
    "radia_windings_qa","radia_material_contract_qa",MaxWindings=4);
verifyEqual(testCase,contract.runtime.winding_count,uint16(1));
verifyTrue(testCase,evalin("base","isa(RadiaWindingBus,'Simulink.Bus')"));
verifyTrue(testCase,evalin("base","isa(RadiaMachineCommandBus,'Simulink.Bus')"));
verifyTrue(testCase,evalin("base","isa(RadiaMachineResponseBus,'Simulink.Bus')"));
end

function testPackagedLibraryContainsWindingBlock(testCase)
repoRoot=fileparts(fileparts(fileparts(mfilename("fullpath"))));
library=fullfile(repoRoot,"matlab","radia_simulink_library.slx");
load_system(library);
cleanup=onCleanup(@() closeIfLoaded("radia_simulink_library")); %#ok<NASGU>
block="radia_simulink_library/Coupling/Winding Dictionary";
verifyEqual(testCase,get_param(block,"BlockType"),'SubSystem');
verifyEqual(testCase,get_param(block,"MaskType"),'Radia Winding Dictionary');
verifyEqual(testCase,get_param(block+"/Compiled Winding Bus","OutDataTypeStr"), ...
    'Bus: RadiaWindingBus');
data=get_param(block,"UserData");
verifyEqual(testCase,data.mechanical_owner,"simulink-or-simscape");
verifyEqual(testCase,data.circuit_backends,"native-mex,ltspice");
end

function contract = materialContract()
names=["air";"phase_a_pos";"phase_a_neg";"phase_b_pos";"phase_b_neg"];
materials=dictionary("copper",{radia.simulink.makeMaterialSpec( ...
    Conductivity_S_per_m=5.8e7)});
regions=dictionary(names,repmat("copper",numel(names),1));
contract=radia.simulink.compileMaterialDictionary(materials, ...
    RegionMaterials=regions,MaxMaterials=4,MaxRegions=8);
end

function cleanupWorkspace
names=["radia_windings_qa","radia_material_contract_qa", ...
    "radia_winding_contract","radia_winding_bus","RadiaWindingBus", ...
    "RadiaMachineCommandBus","RadiaMachineResponseBus", ...
    "radia_machine_command","radia_machine_response"];
for name=names
    if evalin("base","exist('"+name+"','var')"),evalin("base","clear "+name);end
end
end

function closeIfLoaded(name)
if bdIsLoaded(name),close_system(name,0);end
end
