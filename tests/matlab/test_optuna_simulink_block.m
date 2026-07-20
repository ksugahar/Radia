function tests=test_optuna_simulink_block
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
paths=[string(fullfile(root,"matlab")),string(fullfile(root,"tests","matlab","fixtures"))]; addpath(paths(1)); addpath(paths(2)); testCase.TestData.Paths=paths;
end
function teardownOnce(testCase), rmpath(testCase.TestData.Paths(1)); rmpath(testCase.TestData.Paths(2)); end
function testTriggeredBlockRunsStudy(testCase)
m="radia_optuna_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_quadratic",NumTrials=30,SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1','Position',[40 105 100 135]);
add_block('simulink/Sinks/To Workspace',m+"/Best",'VariableName','best_value','SaveFormat','Array','Position',[430 85 520 115]);
add_block('simulink/Sinks/To Workspace',m+"/Status",'VariableName','status_value','SaveFormat','Array','Position',[430 155 520 185]);
add_line(m,'Trigger/1','Optuna Optimization/1'); add_line(m,'Optuna Optimization/1','Best/1'); add_line(m,'Optuna Optimization/3','Status/1');
out=sim(m,'StopTime','0.1','ReturnWorkspaceOutputs','on'); best=out.get('best_value'); status=out.get('status_value');
verifyEqual(testCase,status(end),1); verifyLessThan(testCase,best(end),0.1);
clear cleanup; closeModel(m);
end
function testTriggeredBlockRunsLTspiceStudy(testCase)
m="radia_optuna_ltspice_block_test"; cleanup=onCleanup(@()closeModel(m)); new_system(m);
radia.simulink.buildOptunaBlock(m,ObjectiveFcn="radia_optuna_ltspice_objective",NumTrials=2,SampleTime_s=0.1,Save=false);
add_block('simulink/Sources/Constant',m+"/Trigger",'Value','1'); add_block('simulink/Sinks/To Workspace',m+"/Status",'VariableName','status_value','SaveFormat','Array');
add_line(m,'Trigger/1','Optuna Optimization/1'); add_line(m,'Optuna Optimization/3','Status/1');
out=sim(m,'StopTime','0.1','ReturnWorkspaceOutputs','on'); status=out.get('status_value'); verifyEqual(testCase,status(end),1);
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
function closeModel(m), if bdIsLoaded(m), close_system(m,0); end, end
