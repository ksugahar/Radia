function tests=test_optuna_simulink_block
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
paths=[string(fullfile(root,"matlab")),string(fullfile(root,"tests","matlab","fixtures"))];
entries=string(strsplit(path,pathsep));
removePaths=false(size(paths));
for index=1:numel(paths)
    removePaths(index)=~any(strcmpi(entries,paths(index)));
    if removePaths(index),addpath(paths(index));end
end
testCase.TestData.Paths=paths; testCase.TestData.RemovePaths=removePaths;
end
function teardownOnce(testCase)
for index=1:numel(testCase.TestData.Paths)
    if testCase.TestData.RemovePaths(index),rmpath(testCase.TestData.Paths(index));end
end
end
function testTriggeredBlockRunsStudy(testCase)
m="radia_optuna_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_quadratic",NumTrials=30,SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1','Position',[40 105 100 135]);
add_block('simulink/Sinks/To Workspace',m+"/Best",'VariableName','best_value','SaveFormat','Array','Position',[430 85 520 115]);
add_block('simulink/Sinks/To Workspace',m+"/Status",'VariableName','status_value','SaveFormat','Array','Position',[430 155 520 185]);
add_block('simulink/Sinks/To Workspace',m+"/BestUpdate",'VariableName','best_updated','SaveFormat','Array','Position',[430 215 520 245]);
add_line(m,'Trigger/1','Optuna Optimization/1'); add_line(m,'Optuna Optimization/1','Best/1'); add_line(m,'Optuna Optimization/3','Status/1'); add_line(m,'Optuna Optimization/7','BestUpdate/1');
out=sim(m,'StopTime','3.0','ReturnWorkspaceOutputs','on'); best=out.get('best_value'); status=out.get('status_value'); updates=out.get('best_updated');
verifyEqual(testCase,status(end),1); verifyLessThan(testCase,best(end),0.1);
verifyGreaterThan(testCase,sum(updates),0);
clear cleanup; closeModel(m);
end
function testTriggeredBlockRunsLTspiceStudy(testCase)
m="radia_optuna_ltspice_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_ltspice_objective",NumTrials=2,SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1'); add_block('simulink/Sinks/To Workspace',m+"/Status",'VariableName','status_value','SaveFormat','Array');
add_line(m,'Trigger/1','Optuna Optimization/1'); add_line(m,'Optuna Optimization/3','Status/1');
out=sim(m,'StopTime','0.2','ReturnWorkspaceOutputs','on'); status=out.get('status_value'); verifyEqual(testCase,status(end),1);
clear cleanup; closeModel(m);
end
function testSimulinkOutputsLiveParetoEvolution(testCase)
m="radia_optuna_pareto_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_biobjective",NumTrials=12, ...
    Directions=["minimize","minimize"],SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1');
add_block('simulink/Sinks/To Workspace',m+"/ParetoCount",'VariableName','pareto_count','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/ParetoRevision",'VariableName','pareto_revision','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/BestUpdate",'VariableName','best_updated','SaveFormat','Array');
add_line(m,'Trigger/1','Optuna Optimization/1');
add_line(m,'Optuna Optimization/8','ParetoCount/1');
add_line(m,'Optuna Optimization/11','ParetoRevision/1');
add_line(m,'Optuna Optimization/7','BestUpdate/1');
out=sim(m,'StopTime','1.2','ReturnWorkspaceOutputs','on');
count=out.get('pareto_count'); revision=out.get('pareto_revision'); updates=out.get('best_updated');
verifyGreaterThanOrEqual(testCase,count(end),2);
verifyGreaterThan(testCase,max(revision),1);
verifyGreaterThan(testCase,sum(updates),0);
clear cleanup; closeModel(m);
end
function testMonitorUsesSimulinkScopesAndNoBrowser(testCase)
m="radia_optuna_monitor_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
path=radia.simulink.addOptunaMonitor(m);
verifyEqual(testCase,string(get_param(path+"/Optimization History","BlockType")),"Scope");
verifyEqual(testCase,string(get_param(path+"/Pareto Front","BlockType")),"Record");
contract=get_param(path,"UserData");
verifyFalse(testCase,contract.browser_required);
verifyEqual(testCase,string(contract.visualization),"simulink-scope-xy");
clear cleanup; closeModel(m);
end
function testSheetMetalOptimizationBlockUsesNativeRunner(testCase)
m="radia_sheet_metal_optuna_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
runner=radia.optuna.SheetMetalRunner(@(trial)struct());
[path,monitor]=radia.simulink.buildSheetMetalOptimizationBlock(m, ...
 RunnerVariable="radia_sheet_metal_runner_test",Runner=runner, ...
 NumTrials=7,Directions=["minimize","minimize"],Save=false);
verifyEqual(testCase,string(get_param(path,"FunctionName")),"radia_optuna_sfun");
verifyEqual(testCase,string(get_param(path,"Mask")),"on");
verifyTrue(testCase,contains(string(get_param(path,"Parameters")), ...
    "runner.evaluate(trial)"));
verifyEqual(testCase,string(get_param(path,"num_trials")),"7");
verifyEqual(testCase,string(get_param(path,"sampler_name")),"auto");
ports=get_param(path,"PortHandles");
verifyEqual(testCase,numel(ports.Inport),2);
verifyEqual(testCase,numel(ports.Outport),14);
contract=get_param(path,"UserData");
verifyEqual(testCase,string(contract.domain),"sheet-metal");
verifyEqual(testCase,string(contract.backend),"matlab-native-ngsolve-cubit");
verifyFalse(testCase,contract.browser_required);
verifyEqual(testCase,string(contract.pareto_kernel), ...
    "required-optuna-mex");
verifyFalse(testCase,contract.python_per_trial);
verifyEqual(testCase,string(get_param(monitor,"Mask")),"on");
evalin("base","clear radia_sheet_metal_runner_test");
clear cleanup; closeModel(m);
end
function testObjectiveFunctionHandleParameter(testCase)
m="radia_optuna_handle_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
path=radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_quadratic", ...
 NumTrials=2,SampleTime_s=0.1,Save=false);
set_param(path,"objective_fcn","@(trial)radia_optuna_quadratic(trial)");
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1');
add_block('simulink/Sinks/To Workspace',m+"/Status", ...
 'VariableName','status_value','SaveFormat','Array');
add_line(m,'Trigger/1','Optuna Optimization/1');
add_line(m,'Optuna Optimization/3','Status/1');
out=sim(m,'StopTime','0.2','ReturnWorkspaceOutputs','on');
status=out.get('status_value');
verifyEqual(testCase,status(end),1);
clear cleanup; closeModel(m);
end
function testCAEFailuresAreRecordedAndStudyContinues(testCase)
m="radia_optuna_failure_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
path=radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_mesh_failure", ...
 NumTrials=3,SampleTime_s=0.1,Save=false);
ports=get_param(path,"PortHandles");
verifyEqual(testCase,numel(ports.Inport),2);
verifyEqual(testCase,numel(ports.Outport),14);
contract=get_param(path,"UserData");
verifyEqual(testCase,string(contract.pareto_kernel), ...
    "required-optuna-mex");
verifyEqual(testCase,string(contract.cae_failure_policy),"record-and-continue");
verifyFalse(testCase,contract.python_per_trial);
add_block('simulink/Sources/Constant',m+"/Start",'Value','1');
add_block('simulink/Sources/Constant',m+"/Cancel",'Value','0');
add_block('simulink/Sinks/To Workspace',m+"/Status", ...
 'VariableName','status_value','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/Failed", ...
 'VariableName','failed_value','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/Attempted", ...
 'VariableName','attempted_value','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/FailureCode", ...
 'VariableName','failure_code','SaveFormat','Array');
add_line(m,'Start/1','Optuna Optimization/1');
add_line(m,'Cancel/1','Optuna Optimization/2');
add_line(m,'Optuna Optimization/3','Status/1');
add_line(m,'Optuna Optimization/12','Failed/1');
add_line(m,'Optuna Optimization/13','Attempted/1');
add_line(m,'Optuna Optimization/14','FailureCode/1');
out=sim(m,'StopTime','0.3','ReturnWorkspaceOutputs','on');
status=out.get('status_value'); failed=out.get('failed_value');
attempted=out.get('attempted_value'); code=out.get('failure_code');
verifyEqual(testCase,status(end),1);
verifyEqual(testCase,failed(end),3);
verifyEqual(testCase,attempted(end),3);
verifyEqual(testCase,code(end),3);
clear cleanup; closeModel(m);
end
function testCancelStopsActiveStudy(testCase)
m="radia_optuna_cancel_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_quadratic", ...
 NumTrials=20,SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Start",'Value','1');
add_block('simulink/Sources/Step',m+"/Cancel", ...
 'Time','0.15','Before','0','After','1');
add_block('simulink/Sinks/To Workspace',m+"/Status", ...
 'VariableName','status_value','SaveFormat','Array');
add_block('simulink/Sinks/To Workspace',m+"/Attempted", ...
 'VariableName','attempted_value','SaveFormat','Array');
add_line(m,'Start/1','Optuna Optimization/1');
add_line(m,'Cancel/1','Optuna Optimization/2');
add_line(m,'Optuna Optimization/3','Status/1');
add_line(m,'Optuna Optimization/13','Attempted/1');
out=sim(m,'StopTime','0.5','ReturnWorkspaceOutputs','on');
status=out.get('status_value'); attempted=out.get('attempted_value');
verifyEqual(testCase,status(end),3);
verifyLessThan(testCase,attempted(end),20);
verifyGreaterThan(testCase,attempted(end),0);
clear cleanup; closeModel(m);
end
function testLiveMonitorReceivesTrialProgress(testCase)
monitor=radia.optuna.LiveMonitor(Visible=false); cleanup=onCleanup(@()delete(monitor));
study=radia.optuna.createStudy(AutoSave=false,ProgressFcn=@monitor.update);
study.optimize(@radia_optuna_quadratic,3);
verifyEqual(testCase,monitor.UpdateCount,6);
verifyTrue(testCase,isgraphics(monitor.Figure));
clear cleanup; delete(monitor);
end
function testLiveMonitorDisplaysParetoProgress(testCase)
monitor=radia.optuna.LiveMonitor(Visible=false); cleanup=onCleanup(@()delete(monitor));
study=radia.optuna.createStudy(directions=["minimize","minimize"], ...
    AutoSave=false,ProgressFcn=@monitor.update);
for values=[0 2;1 1;2 0]'
    trial=study.ask(); study.tell(trial,values');
end
verifyEqual(testCase,monitor.UpdateCount,6);
verifyEqual(testCase,height(study.paretoFront()),3);
clear cleanup; delete(monitor);
end
function testBuilderExposesLargeNativeSamplerChoices(testCase)
m="radia_optuna_large_sampler_choices";
cleanup=onCleanup(@()closeModel(m));
new_system(m);
for choice=["gp","nsgaiii","bruteforce","qmc"]
    path=radia.simulink.buildOptunaBlock(m, ...
        ObjectiveFcn="radia_optuna_quadratic",NumTrials=2, ...
        Sampler=choice,Save=false);
    verifyEqual(testCase,string(get_param(path,"sampler_name")),choice);
    parameter=Simulink.Mask.get(path).getParameter("sampler_name");
    options=string(parameter.TypeOptions);
    verifyTrue(testCase,all(ismember( ...
        ["gp","nsgaiii","bruteforce","qmc"],options)));
    delete_block(path);
end
clear cleanup; closeModel(m);
end
function closeModel(m), if bdIsLoaded(m), close_system(m,0); end, end
