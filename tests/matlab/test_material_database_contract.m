function tests = test_material_database_contract
%TEST_MATERIAL_DATABASE_CONTRACT QA for spatial material data shared by S-Functions.
tests = functiontests(localfunctions);
end

function testValidLocalEvaluation(testCase)
db = fixtureDatabase();
state = radia.material.makeFieldState([293.15;343.15],[0;50],MeshId="mesh-1");
local = radia.material.evaluateLocal(db,state);
verifyEqual(testCase,local.coordinate_system,"workpiece");
verifySize(testCase,local.B_T,[2 1]);
verifyGreaterThan(testCase,min(local.dBdH_T_per_Apm),0);
verifyGreaterThan(testCase,min(local.conductivity_S_per_m),0);
end

function testFieldStateRejectsDifferentSizes(testCase)
verifyError(testCase,@() radia.material.makeFieldState([293.15;300],[1]), ...
    "radia:material:FieldStateSize");
end

function testDatabaseRejectsNonpositiveDifferentialSlope(testCase)
db = fixtureDatabase(); db.bh_dBdH(1,1) = 0;
verifyError(testCase,@() radia.material.validateDatabase(db),"radia:material:BH");
end

function testLibraryRegistersMaterialDatabaseBlock(testCase)
load_system("simulink");
root = fullfile(tempdir,"radia_material_block_test");
if isfolder(root), rmdir(root,"s"); end
mkdir(root);
addpath(fileparts(fileparts(fileparts(mfilename("fullpath")))),'-begin');
addpath("C:/temp/radia_ih_native_mex2",'-begin');
path = radia.simulink.buildLibrary(OutputDirectory=root);
cleanup = onCleanup(@() closeIfLoaded("radia_simulink_library")); %#ok<NASGU>
load_system(path);
block = "radia_simulink_library/Material Models/Material Database";
verifyEqual(testCase,get_param(block,"BlockType"),'SubSystem');
verifyNotEmpty(testCase,Simulink.Mask.get(block));
data = get_param(block,"UserData");
verifyEqual(testCase,data.schema,"radia.material.database.v1");
end

function db = fixtureDatabase()
db = struct("schema","radia.material.database.v1", ...
    "material_id","qa-steel", "coordinate_system","workpiece", ...
    "conductivity_S_per_m",[1.0e6 0.9e6;0.8e6 0.7e6], ...
    "bh_temperature_K",[293.15;393.15], ...
    "bh_H_A_per_m",[0;100], ...
    "bh_B_T",[0 0.1;0 0.08], ...
    "bh_dBdH",[1e-3 8e-4;9e-4 7e-4]);
end

function closeIfLoaded(name)
if bdIsLoaded(name), close_system(name,0); end
end
