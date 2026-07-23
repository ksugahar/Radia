function tests = test_bh_sfun
%TEST_BH_SFUN Unit and library-contract tests for the temperature BH block.
tests = functiontests(localfunctions);
end

function testFormulaBlock(testCase)
cfg = struct("mode","formula", "reference_temperature_K",293.15, ...
    "mu_r_ref",100, "mu_r_temperature_slope",-0.1);
out = runBHModel(cfg,293.15,2.0);
verifyEqual(testCase,out.B,4*pi*1e-7*100*2,"AbsTol",1e-12);
verifyEqual(testCase,out.dBdH,4*pi*1e-7*100,"AbsTol",1e-12);
end

function testLUTBlock(testCase)
cfg = struct("mode","lut", "temperature_K",[293.15;393.15], ...
    "H_A_per_m",[0;100], "B_T",[0 0.1;0 0.2]);
out = runBHModel(cfg,343.15,50);
verifyEqual(testCase,out.B,0.075,"AbsTol",1e-12);
verifyGreaterThan(testCase,out.dBdH,0);
end

function testLibraryRegistersBH(testCase)
load_system("simulink");
root = fullfile("C:\temp","radia_bh_test_library");
if isfolder(root), rmdir(root,"s"); end
mkdir(root);
addpath(fileparts(fileparts(fileparts(mfilename("fullpath")))),'-begin');
path = radia.simulink.buildLibrary(OutputDirectory=root);
cleanup = onCleanup(@() cleanupLibrary(root));
load_system(path);
block = "radia_simulink_library/Material Models/Temperature-Dependent BH";
verifyTrue(testCase,bdIsLoaded("radia_simulink_library"));
verifyEqual(testCase,get_param(block,"FunctionName"),'radia_bh_sfun');
verifyNotEmpty(testCase,Simulink.Mask.get(block));
end

function out = runBHModel(cfg,T,H)
model = "radia_bh_unit_" + erase(string(java.util.UUID.randomUUID),"-");
assignin("base","radia_bh_test_config",cfg);
new_system(model);
cleanup = onCleanup(@() closeIfLoaded(model));
add_block("simulink/Sources/Constant",model+"/T","Value",num2str(T),"Position",[30 40 80 70]);
add_block("simulink/Sources/Constant",model+"/H","Value",num2str(H),"Position",[30 120 80 150]);
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function",model+"/BH", ...
    "FunctionName","radia_bh_sfun","Parameters","radia_bh_test_config", ...
    "Position",[140 55 280 135]);
add_block("simulink/Sinks/To Workspace",model+"/B","VariableName","bh_B","SaveFormat","Array", ...
    "Position",[340 45 430 75]);
add_block("simulink/Sinks/To Workspace",model+"/dBdH","VariableName","bh_dBdH","SaveFormat","Array", ...
    "Position",[340 105 430 135]);
add_line(model,"T/1","BH/1"); add_line(model,"H/1","BH/2");
add_line(model,"BH/1","B/1"); add_line(model,"BH/2","dBdH/1");
set_param(model,"StopTime","1");
simOut = sim(model,"ReturnWorkspaceOutputs","on");
out.B = simOut.get("bh_B"); out.B = out.B(end);
out.dBdH = simOut.get("bh_dBdH"); out.dBdH = out.dBdH(end);
end

function closeIfLoaded(name)
if bdIsLoaded(name), close_system(name,0); end
end

function cleanupLibrary(root)
closeIfLoaded("radia_simulink_library");
if isfolder(root), rmdir(root,"s"); end
end
