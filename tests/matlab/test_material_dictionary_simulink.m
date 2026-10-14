function tests = test_material_dictionary_simulink
%TEST_MATERIAL_DICTIONARY_SIMULINK Fixed-width material and .vol contracts.
tests = functiontests(localfunctions);
end

function testCompilesDictionaryAgainstVol(testCase)
root = freshRoot();
cleanup = onCleanup(@() removeRoot(root)); %#ok<NASGU>
mesh = writeVol(root,["air_domain","stator","phase_a","rotor_pm"]);

air = radia.simulink.makeMaterialSpec();
steel = radia.simulink.makeMaterialSpec( ...
    BH_B_T=[0;1.0;1.5],BH_H_A_per_m=[0;250;1200]);
copper = radia.simulink.makeMaterialSpec( ...
    Conductivity_S_per_m=5.8e7,Density_kg_per_m3=8960, ...
    SpecificHeat_J_per_kgK=385,ThermalConductivity_W_per_mK=401);
magnet = radia.simulink.makeMaterialSpec(MuR=1.05,Remanence_T=[0 0 1.2]);
materials = dictionary(["air","steel","copper","magnet"], ...
    {air,steel,copper,magnet});
regions = dictionary(["air_domain","stator","phase_a","rotor_pm"], ...
    ["air","steel","copper","magnet"]);

contract = radia.simulink.compileMaterialDictionary(materials, ...
    RegionMaterials=regions,MeshFile=mesh,MaxMaterials=8,MaxRegions=8, ...
    MaxBHPoints=8,MaxHysteresisParameters=4);
verifyEqual(testCase,contract.schema,"radia.simulink.material-dictionary.v1");
verifyEqual(testCase,contract.material_names,["air";"copper";"magnet";"steel"]);
verifyEqual(testCase,contract.runtime.material_count,uint16(4));
verifyEqual(testCase,contract.runtime.region_count,uint16(4));
verifyEqual(testCase,contract.runtime.region_material_index(1:4),uint16([1;4;2;3]));
verifyEqual(testCase,contract.runtime.bh_count(4),uint16(3));
verifyNotEqual(testCase,bitand(contract.runtime.kind_flags(2),uint16(4)),uint16(0));
verifyNotEqual(testCase,bitand(contract.runtime.kind_flags(3),uint16(8)),uint16(0));
verifyFalse(testCase,contract.runtime_policy.dictionary_lookup_per_step);
verifyFalse(testCase,contract.runtime_policy.python_per_step);
verifyEqual(testCase,strlength(contract.mesh.mesh_sha256),64);

runtimeFields = fieldnames(contract.runtime);
verifyTrue(testCase,all(cellfun(@(name) isnumeric(contract.runtime.(name)) || ...
    islogical(contract.runtime.(name)),runtimeFields)));
bus = radia.simulink.makeMaterialBusObject(contract.runtime,AssignToBase=false);
verifyEqual(testCase,numel(bus.Elements),numel(runtimeFields));
end

function testSameNameMappingNeedsNoRegionDictionary(testCase)
root = freshRoot();
cleanup = onCleanup(@() removeRoot(root)); %#ok<NASGU>
mesh = writeVol(root,["air","steel"]);
materials = dictionary(["air","steel"], ...
    {radia.simulink.makeMaterialSpec(), ...
     radia.simulink.makeMaterialSpec(MuR=1000)});
contract = radia.simulink.compileMaterialDictionary(materials,MeshFile=mesh);
verifyEqual(testCase,contract.region_material_names,["air";"steel"]);
verifyEqual(testCase,contract.runtime.region_material_index(1:2),uint16([1;2]));
end

function testRejectsIncompleteRegionCoverage(testCase)
root = freshRoot();
cleanup = onCleanup(@() removeRoot(root)); %#ok<NASGU>
mesh = writeVol(root,["air_domain","core"]);
materials = dictionary(["air","steel"], ...
    {radia.simulink.makeMaterialSpec(),radia.simulink.makeMaterialSpec(MuR=1000)});
regions = dictionary("air_domain","air");
verifyError(testCase,@() radia.simulink.compileMaterialDictionary(materials, ...
    RegionMaterials=regions,MeshFile=mesh),"radia:simulink:RegionCoverage");
end

function testRejectsUnknownRegionMaterial(testCase)
root = freshRoot();
cleanup = onCleanup(@() removeRoot(root)); %#ok<NASGU>
mesh = writeVol(root,"core");
materials = dictionary("steel",{radia.simulink.makeMaterialSpec(MuR=1000)});
regions = dictionary("core","missing");
verifyError(testCase,@() radia.simulink.compileMaterialDictionary(materials, ...
    RegionMaterials=regions,MeshFile=mesh),"radia:simulink:UnknownRegionMaterial");
end

function testRejectsCapacityAndStaleMesh(testCase)
root = freshRoot();
cleanup = onCleanup(@() removeRoot(root)); %#ok<NASGU>
mesh = writeVol(root,["air","steel"]);
materials = dictionary(["air","steel"], ...
    {radia.simulink.makeMaterialSpec(),radia.simulink.makeMaterialSpec(MuR=1000)});
verifyError(testCase,@() radia.simulink.compileMaterialDictionary(materials, ...
    MeshFile=mesh,MaxMaterials=1),"radia:simulink:MaterialCapacity");
verifyError(testCase,@() radia.simulink.compileMaterialDictionary(materials, ...
    MeshFile=mesh,ExpectedMeshSHA256=string(repmat('0',1,64))), ...
    "radia:simulink:MeshDigest");
end

function testRejectsBadConstitutiveData(testCase)
verifyError(testCase,@() radia.simulink.makeMaterialSpec( ...
    BH_B_T=[0;1;0.9],BH_H_A_per_m=[0;100;200]), ...
    "radia:simulink:MaterialBHMonotonic");
verifyError(testCase,@() radia.simulink.makeMaterialSpec( ...
    Density_kg_per_m3=7800),"radia:simulink:MaterialThermalTuple");
verifyError(testCase,@() radia.simulink.makeMaterialSpec( ...
    HysteresisModel="energy"),"radia:simulink:MaterialHysteresisParameters");
end

function testRegistersWorkspaceDictionaryAndBus(testCase)
root = freshRoot();
cleanup = onCleanup(@() cleanupWorkspace(root)); %#ok<NASGU>
mesh = writeVol(root,["air","coil"]);
materials = dictionary(["air","copper"], ...
    {radia.simulink.makeMaterialSpec(), ...
     radia.simulink.makeMaterialSpec(Conductivity_S_per_m=5.8e7)});
regions = dictionary(["air","coil"],["air","copper"]);
assignin("base","radia_materials_qa",materials);
assignin("base","radia_regions_qa",regions);
contract = radia.simulink.registerMaterialDictionary( ...
    "radia_materials_qa","radia_regions_qa",mesh,ValidateMesh=false, ...
    MaxMaterials=4,MaxRegions=4,MaxBHPoints=8,MaxHysteresisParameters=4);
verifyEqual(testCase,contract.runtime.region_material_index(1:2),uint16([1;2]));
verifyTrue(testCase,evalin("base","isa(RadiaMaterialBus,'Simulink.Bus')"));
verifyTrue(testCase,evalin("base","isstruct(radia_material_bus)"));
end

function testPackagedLibraryContainsTypedMaterialBlock(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
library = fullfile(repoRoot,"matlab","radia_simulink_library.slx");
load_system(library);
cleanup = onCleanup(@() closeIfLoaded("radia_simulink_library")); %#ok<NASGU>
block = "radia_simulink_library/Material Models/Material Dictionary";
verifyEqual(testCase,get_param(block,"BlockType"),'SubSystem');
verifyEqual(testCase,get_param(block,"MaskType"),'Radia Material Dictionary');
verifyEqual(testCase,get_param(block+"/Compiled Material Bus","OutDataTypeStr"), ...
    'Bus: RadiaMaterialBus');
data = get_param(block,"UserData");
verifyEqual(testCase,data.schema,"radia.simulink.material-dictionary.v1");
end

function testLibraryContainsTypedMaterialBlock(testCase)
assumeTrue(testCase,exist("radia_mex","file") == 3 && ...
    exist("radia_ih_eddy_sfun","file") == 2, ...
    "The IH MEX ABI and Level-2 MATLAB wrappers are required to rebuild the library.");
root = freshRoot();
cleanup = onCleanup(@() cleanupLibrary(root)); %#ok<NASGU>
path = radia.simulink.buildLibrary(OutputDirectory=root);
load_system(path);
block = "radia_simulink_library/Material Models/Material Dictionary";
verifyEqual(testCase,get_param(block,"BlockType"),'SubSystem');
verifyNotEmpty(testCase,Simulink.Mask.get(block));
data = get_param(block,"UserData");
verifyEqual(testCase,data.schema,"radia.simulink.material-dictionary.v1");
verifyEqual(testCase,data.runtime_bus,"RadiaMaterialBus");
verifyFalse(testCase,data.dictionary_lookup_per_step);
verifyEqual(testCase,get_param(block+"/Compiled Material Bus","OutDataTypeStr"), ...
    'Bus: RadiaMaterialBus');
end

function root = freshRoot()
root = string(tempname("C:\temp"));
mkdir(root);
end

function mesh = writeVol(root,names)
names = string(names(:));
mesh = fullfile(root,"materials.vol");
file = fopen(mesh,"wt");
cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
fprintf(file,"mesh3d\ndimension\n3\nmaterials\n%d\n",numel(names));
for k = 1:numel(names), fprintf(file,"%d %s\n",k,names(k)); end
fprintf(file,"volumeelements\n0\nsurfaceelements\n0\nedgesegmentsgi2\n0\npoints\n0\n");
end

function cleanupLibrary(root)
if bdIsLoaded("radia_simulink_library"),close_system("radia_simulink_library",0);end
removeRoot(root);
end

function cleanupWorkspace(root)
variables = ["radia_materials_qa","radia_regions_qa","radia_material_contract", ...
    "radia_material_bus","RadiaMaterialBus"];
for name = variables
    if evalin("base","exist('"+name+"','var')"),evalin("base","clear "+name);end
end
removeRoot(root);
end

function closeIfLoaded(name)
if bdIsLoaded(name),close_system(name,0);end
end

function removeRoot(root)
if isfolder(root),rmdir(root,"s");end
end
