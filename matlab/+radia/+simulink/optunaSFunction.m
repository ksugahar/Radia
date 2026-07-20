function optunaSFunction(block)
%OPTUNASFUNCTION Level-2 block that starts a MATLAB Optuna study on trigger.
setup(block);
end
function setup(b)
b.NumDialogPrms=6; b.DialogPrmsTunable=repmat({'Nontunable'},1,6);
b.NumInputPorts=1; b.NumOutputPorts=6; b.SetPreCompInpPortInfoToDynamic; b.SetPreCompOutPortInfoToDynamic;
b.InputPort(1).Dimensions=1; b.InputPort(1).DirectFeedthrough=true;
for k=1:6, b.OutputPort(k).Dimensions=1; b.OutputPort(k).DatatypeID=0; end
b.SampleTimes=[b.DialogPrm(5).Data 0];
b.RegBlockMethod('PostPropagationSetup',@postSetup); b.RegBlockMethod('Start',@start); b.RegBlockMethod('Outputs',@outputs);
end
function postSetup(b)
b.NumDworks=7; names={'previous_trigger','best_value','best_trial','status','completed_trials','last_value','elapsed_s'};
for k=1:7
 b.Dwork(k).Name=names{k}; b.Dwork(k).Dimensions=1; b.Dwork(k).DatatypeID=0; b.Dwork(k).Complexity='Real'; b.Dwork(k).UsedAsDiscState=true;
end
end
function start(b)
for k=1:7, b.Dwork(k).Data=0; end
b.Dwork(2).Data=NaN; b.Dwork(3).Data=NaN; b.Dwork(6).Data=NaN;
end
function outputs(b)
trigger=double(b.InputPort(1).Data);
if trigger>0 && b.Dwork(1).Data<=0
 objective=str2func(char(string(b.DialogPrm(1).Data))); nTrials=double(b.DialogPrm(2).Data);
 direction=string(b.DialogPrm(3).Data); storage=string(b.DialogPrm(4).Data);
 liveVisualization=logical(b.DialogPrm(6).Data); started=tic;
 try
  monitor=[];
  if liveVisualization, monitor=radia.optuna.LiveMonitor(); end
  progress=[]; if ~isempty(monitor), progress=@monitor.update; end
  study=radia.optuna.createStudy(direction=direction,StoragePath=storage,AutoSave=strlength(storage)>0,ProgressFcn=progress);
  study.optimize(objective,nTrials); best=study.bestTrial();
  b.Dwork(2).Data=best.Value(1); b.Dwork(3).Data=best.TrialNumber(1); b.Dwork(4).Data=1;
  complete=study.TrialTable.State=="COMPLETE"; b.Dwork(5).Data=sum(complete);
  if any(complete), b.Dwork(6).Data=study.TrialTable.Value(find(complete,1,"last")); end
  b.Dwork(7).Data=toc(started);
 catch exception
  b.Dwork(4).Data=-1; warning('radia:simulink:OptunaFailed','Optuna block failed: %s',exception.message);
 end
end
b.Dwork(1).Data=trigger;
for k=1:6, b.OutputPort(k).Data=b.Dwork(k+1).Data; end
end
