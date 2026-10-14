function tests = test_electromagnet_topology_model
%TEST_ELECTROMAGNET_TOPOLOGY_MODEL Verify the VIM density-adjoint model.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDirectory = fileparts(mfilename("fullpath"));
repositoryRoot = fileparts(fileparts(testDirectory));
matlabDirectory = fullfile(repositoryRoot,"matlab");
entries = string(strsplit(path,pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
testCase.TestData.FileGenConfig = Simulink.fileGenControl("getConfig");
testCase.TestData.FileGenRoot = string(tempname("C:\temp"));
Simulink.fileGenControl("set", ...
    CacheFolder=fullfile(testCase.TestData.FileGenRoot,"cache"), ...
    CodeGenFolder=fullfile(testCase.TestData.FileGenRoot,"codegen"), ...
    createDir=true);
end

function teardownOnce(testCase)
Simulink.fileGenControl("setConfig", ...
    config=testCase.TestData.FileGenConfig);
if isfolder(testCase.TestData.FileGenRoot)
    rmdir(testCase.TestData.FileGenRoot,"s");
end
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testReferenceRunnerHasCheckedAdjointGradient(testCase)
runner = radia.topopt.makeElectromagnetReferenceRunner();
verifyEqual(testCase,string(runner.Metadata.domain), ...
    "electromagnet-topology");
verifyEqual(testCase,string(runner.Metadata.gradient), ...
    "one-state-one-adjoint");
diagnostic = radia.topopt.checkAdjointGradient( ...
    runner.InitialDesign,runner.EvaluateFcn,RelativeStep=2e-6);
verifyTrue(testCase,diagnostic.passed);
verifyLessThan(testCase,max(diagnostic.history.ObjectiveError),2e-7);
verifyLessThan(testCase,max(diagnostic.history.ConstraintError),1e-9);
end

function testReferenceRunnerOptimizesDensityAndVolume(testCase)
assumeTrue(testCase,hasOptimizationToolbox(), ...
    "Optimization Toolbox is required for MMA/SQP.");
runner = radia.topopt.makeElectromagnetReferenceRunner();
initial = runner.EvaluateFcn(runner.InitialDesign);
result = runner.run();
verifyTrue(testCase,result.converged,result.output.message);
verifyLessThan(testCase,result.objective,initial.objective);
verifyLessThanOrEqual(testCase,max(result.constraints),1e-7);
verifyGreaterThanOrEqual(testCase,min(result.design),0);
verifyLessThanOrEqual(testCase,max(result.design),1);
verifyEqual(testCase,result.evaluation.payload.volume_fraction, ...
    0.5,"AbsTol",5e-4);
end

function testReferenceRunnerSupportsSQP(testCase)
assumeTrue(testCase,hasOptimizationToolbox(), ...
    "Optimization Toolbox is required for MMA/SQP.");
runner = radia.topopt.makeElectromagnetReferenceRunner(Solver="sqp");
result = runner.run();
verifyTrue(testCase,result.converged,result.output.message);
verifyEqual(testCase,result.solver,"sqp");
verifyLessThan(testCase,result.objective,1e-7);
verifyLessThanOrEqual(testCase,max(result.constraints),1e-8);
verifyEqual(testCase,result.evaluation.payload.volume_fraction, ...
    0.5,"AbsTol",1e-8);
end

function testBuilderCreatesRunnableModel(testCase)
outputDirectory = string(tempname("C:\temp"));
mkdir(outputDirectory);
modelName = "radia_electromagnet_test_" + ...
    erase(string(java.util.UUID.randomUUID),"-");
cleanup = onCleanup(@() closeIfLoaded(modelName));
modelPath = radia.simulink.buildElectromagnetOptimizationModel( ...
    ModelName=modelName,OutputDirectory=outputDirectory, ...
    SampleTime_s=0.1,HistoryCapacity=96,Open=false);
verifyTrue(testCase,isfile(modelPath));
load_system(modelPath);
blockPath = modelName + "/Electromagnet Topology Optimization";
verifyEqual(testCase,string(get_param(blockPath,"Mask")),"on");
verifyEqual(testCase,string(get_param( ...
    blockPath + "/Density Topology Optimization","FunctionName")), ...
    "radia_electromagnet_topology_sfun");
contract = get_param(blockPath,"UserData");
verifyEqual(testCase,string(contract.domain), ...
    "electromagnet-topology-optimization");
verifyEqual(testCase,contract.state_solves_per_gradient,1);
verifyEqual(testCase,contract.adjoint_solves_per_gradient,1);
verifyFalse(testCase,contract.finite_difference_fallback);
verifyFalse(testCase,contract.python_per_optimization_step);
set_param(modelName,"SimulationCommand","update");
simulation = sim(modelName,"ReturnWorkspaceOutputs","on");
outputs = simulation.get("yout");
status = outputs.getElement(5).Values.Data;
objective = outputs.getElement(4).Values.Data;
density = outputs.getElement(7).Values.Data;
volumeFraction = outputs.getElement(8).Values.Data;
if hasOptimizationToolbox()
    verifyEqual(testCase,status(end),2);
    verifyTrue(testCase,isfinite(objective(end)));
    verifySize(testCase,density,[3,4]);
    verifyLessThanOrEqual(testCase,volumeFraction(end),0.5+1e-7);
else
    verifyEqual(testCase,status(end),-1);
end
clear cleanup
closeIfLoaded(modelName);
end

function testTrackedModelLoadsAndUpdates(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
modelPath = fullfile(root,"matlab","radia_electromagnet.slx");
verifyTrue(testCase,isfile(modelPath));
load_system(modelPath);
cleanup = onCleanup(@() closeIfLoaded("radia_electromagnet"));
set_param("radia_electromagnet","SimulationCommand","update");
blockPath = "radia_electromagnet/Electromagnet Topology Optimization";
verifyEqual(testCase,string(get_param( ...
    blockPath + "/Density Topology Optimization","FunctionName")), ...
    "radia_electromagnet_topology_sfun");
verifyEqual(testCase,string(get_param(blockPath,"Tag")), ...
    "RadiaElectromagnetTopology");
clear cleanup
closeIfLoaded("radia_electromagnet");
end

function available = hasOptimizationToolbox()
available = ~isempty(ver("optim")) && ...
    license("test","Optimization_Toolbox") && ...
    exist("fmincon","file") == 2;
end

function closeIfLoaded(name)
if bdIsLoaded(name)
    close_system(name,0);
end
end
