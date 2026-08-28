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
fixture=fullfile(root,"tests","matlab","fixtures","optuna49_oracle.json");
testCase.TestData.Oracle=jsondecode(fileread(fixture));
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testEnqueueTrialMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.core_api.enqueue;
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(7), ...
    AutoSave=false);
params=struct("x",0.25,"mesh",5,"mode","B");
study.enqueue_trial(params,UserAttrs=struct("source","baseline"));
study.enqueue_trial(params,SkipIfExists=true);
verifyEqual(testCase,height(study.TrialTable), ...
    double(expected.waiting_count_after_skip));
trial=study.ask();
actual=struct( ...
    "number",trial.Number, ...
    "x",trial.suggestFloat("x",0,1), ...
    "mesh",trial.suggest_int("mesh",1,9,Step=2), ...
    "mode",trial.suggestCategorical("mode",["A","B"]), ...
    "source",string(trial.UserAttrs.source));
study.tell(trial,1);
verifyEqual(testCase,actual.number,double(expected.values.number));
verifyEqual(testCase,actual.x,double(expected.values.x),AbsTol=0);
verifyEqual(testCase,actual.mesh,double(expected.values.mesh));
verifyEqual(testCase,string(actual.mode),string(expected.values.mode));
verifyEqual(testCase,actual.source,string(expected.values.source));
verifyEqual(testCase,study.TrialTable.State,string(expected.final_state));
end

function testInvalidEnqueuedValueFallsBackLikeUpstream(testCase)
expected=double(testCase.TestData.Oracle.core_api. ...
    invalid_enqueued_fallback);
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(9), ...
    AutoSave=false);
study.enqueueTrial(struct("x",2));
trial=study.ask();
warning("off","radia:optuna:FixedParameter");
cleanup=onCleanup(@()warning("on","radia:optuna:FixedParameter"));
actual=trial.suggestFloat("x",0,1);
verifyEqual(testCase,actual,expected,AbsTol=0);
clear cleanup
end

function testAskFixedDistributionsMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.core_api.ask_fixed_distributions;
distributions=struct( ...
    "x",radia.optuna.FloatDistribution(-1,1), ...
    "mesh",radia.optuna.IntDistribution(1,9,Step=2), ...
    "mode",radia.optuna.CategoricalDistribution(["A","B"]));
study=radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(19),AutoSave=false);
trial=study.ask(distributions);
verifyEqual(testCase,trial.Params.x,double(expected.x),AbsTol=0);
verifyEqual(testCase,trial.Params.mesh,double(expected.mesh));
verifyEqual(testCase,string(trial.Params.mode),string(expected.mode));
verifyEqual(testCase,trial.suggest_float("x",-1,1),double(expected.x), ...
    AbsTol=0);
verifyEqual(testCase,trial.suggest_int("mesh",1,9,Step=2), ...
    double(expected.mesh));
verifyEqual(testCase,string(trial.suggest_categorical( ...
    "mode",["A","B"])),string(expected.mode));
end
