function tests=test_optuna_core_parity
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
testCase.TestData.RemovePath=~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemovePath, addpath(matlabDirectory); end
testCase.TestData.MatlabDirectory=matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testEnqueueTrialIsFifoAndOverridesSampler(testCase)
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(7), ...
    AutoSave=false);
params=struct("x",0.25,"mesh",5,"mode","B");
study.enqueue_trial(params,UserAttrs=struct("source","baseline"));
study.enqueue_trial(params,SkipIfExists=true);
verifyEqual(testCase,height(study.TrialTable),1);
verifyEqual(testCase,study.TrialTable.State,"WAITING");
trial=study.ask();
verifyEqual(testCase,trial.Number,0);
verifyEqual(testCase,trial.suggestFloat("x",0,1),0.25);
verifyEqual(testCase,trial.suggest_int("mesh",1,9,Step=2),5);
verifyEqual(testCase,trial.suggestCategorical("mode",["A","B"]),"B");
verifyEqual(testCase,string(trial.UserAttrs.source),"baseline");
study.tell(trial,1);
verifyEqual(testCase,study.TrialTable.State,"COMPLETE");
end

function testInvalidEnqueuedValueFallsBackToSampler(testCase)
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(9), ...
    AutoSave=false);
study.enqueueTrial(struct("x",2));
trial=study.ask();
lastwarn("");
x=trial.suggestFloat("x",0,1);
[~,warningId]=lastwarn;
verifyEqual(testCase,string(warningId),"radia:optuna:FixedParameter");
verifyGreaterThanOrEqual(testCase,x,0);
verifyLessThanOrEqual(testCase,x,1);
end

function testOptimizeCatchCallbackAndStop(testCase)
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(11), ...
    AutoSave=false);
study.setSystemAttr("callback_count",0);
results=study.optimize(@(trial)controlledObjective(trial,study),10, ...
    Catch="radia:test:Expected",Callbacks=@countCallback);
verifyEqual(testCase,height(results),2);
verifyEqual(testCase,results.State,["FAIL";"COMPLETE"]);
verifyEqual(testCase,study.SystemAttrs.callback_count,2);
verifyEqual(testCase,height(study.ObjectiveTable),1);
end

function testZeroTimeoutStartsNoTrial(testCase)
study=radia.optuna.Study(AutoSave=false);
study.optimize(@(~)0,100,Timeout=0);
verifyEmpty(testCase,study.TrialTable);
verifyError(testCase,@()study.stop(), ...
    "radia:optuna:StopOutsideOptimize");
end

function testMetricSystemAttrsAndFrozenTrial(testCase)
study=radia.optuna.Study(Directions=["minimize","maximize"], ...
    AutoSave=false);
study.set_metric_names(["loss","strength"]);
study.set_system_attr("owner","cae");
trial=study.ask();
trial.set_user_attr("tag","reference");
trial.set_system_attr("worker",3);
trial.suggestFloat("x",0,1);
study.tell(trial,[1,2]);
frozen=study.freezeTrial(trial.Number);
verifyEqual(testCase,study.metric_names(),["loss","strength"]);
verifyEqual(testCase,study.system_attrs().owner,"cae");
verifyClass(testCase,frozen,"radia.optuna.FrozenTrial");
verifyEqual(testCase,frozen.Values,[1,2]);
verifyEqual(testCase,string(frozen.UserAttrs.tag),"reference");
verifyEqual(testCase,frozen.SystemAttrs.worker,3);
end

function testAddTrialImportsDistributionsAndConstraints(testCase)
params=struct("x",0.5,"mesh",5,"mode","B");
distributions=struct( ...
    "x",radia.optuna.FloatDistribution(0,1), ...
    "mesh",radia.optuna.IntDistribution(1,9,Step=2), ...
    "mode",radia.optuna.CategoricalDistribution(["A","B"]));
frozen=radia.optuna.createTrial(Values=[1,2],Params=params, ...
    Distributions=distributions,Constraints=zeros(1,0), ...
    ConstraintPresent=true,UserAttrs=struct("source","archive"));
study=radia.optuna.Study(Directions=["minimize","minimize"], ...
    AutoSave=false);
study.add_trial(frozen);
verifyEqual(testCase,study.TrialTable.State,"COMPLETE");
verifyEqual(testCase,height(study.ParamTable),3);
verifyEqual(testCase,height(study.ObjectiveTable),2);
[present,constraints]=study.constraintRecord(0);
verifyTrue(testCase,present);
verifyEmpty(testCase,constraints);
verifyEqual(testCase,string(study.freezeTrial(0).UserAttrs.source),"archive");
end

function value=controlledObjective(trial,study)
trial.suggestFloat("x",0,1);
if trial.Number==0
    throw(MException("radia:test:Expected","expected objective failure"));
end
study.stop();
value=trial.Number;
end

function countCallback(study,~)
count=study.SystemAttrs.callback_count;
study.setSystemAttr("callback_count",count+1);
end
