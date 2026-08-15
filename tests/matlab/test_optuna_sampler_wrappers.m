function tests=test_optuna_sampler_wrappers
tests=functiontests(localfunctions);
end

function testPartialFixedDelegatesOnlyUnfixedParameters(testCase)
base=radia.optuna.RandomSampler(41);
sampler=radia.optuna.PartialFixedSampler(struct("x",0.25),base);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
trial=study.ask();
x=trial.suggestFloat("x",0,1);
y=trial.suggestFloat("y",0,1);
verifyEqual(testCase,x,0.25);
verifyGreaterThanOrEqual(testCase,y,0);
verifyLessThanOrEqual(testCase,y,1);
study.tell(trial,x+y);
verifyEqual(testCase,study.bestParams().x,0.25);
end

function testPartialFixedKeepsOptunaOutOfRangeWarningBehavior(testCase)
sampler=radia.optuna.PartialFixedSampler( ...
    struct("x",2),radia.optuna.RandomSampler(2));
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
trial=study.ask();
verifyWarning(testCase,@()trial.suggestFloat("x",0,1), ...
    "radia:optuna:FixedParameter");
verifyEqual(testCase,trial.Params.x,2);
end

function testGridSamplerCoversCartesianProductAndStops(testCase)
space=struct("x",[1 2],"mode",{{"a","b"}});
sampler=radia.optuna.GridSampler(space,Seed=3);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
result=study.optimize(@gridObjective,10);
verifyEqual(testCase,height(result),4);
verifyTrue(testCase,sampler.isExhausted(study));
keys=strings(1,height(result));
for row=1:height(result)
    params=result.Params{row};
    keys(row)=string(params.x)+":"+string(params.mode);
end
verifyEqual(testCase,numel(unique(keys)),4);
end

function testSobolSequenceMatchesOptunaUnscrambledPrefix(testCase)
sampler=radia.optuna.QMCSampler(QMCType="sobol",Seed=7);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
addNumericTrial(study); % Initial trial is sampled independently.
points=collectNumericPoints(study,5);
expected=[0 0 0;0.5 0.5 0.5;0.75 0.25 0.25; ...
    0.25 0.75 0.75;0.375 0.375 0.625];
verifyEqual(testCase,points,expected,"AbsTol",0);
end

function testHaltonSequenceMatchesOptunaUnscrambledPrefix(testCase)
sampler=radia.optuna.QMCSampler(QMCType="halton",Seed=7);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
addNumericTrial(study);
points=collectNumericPoints(study,5);
expected=[0 0 0;0.5 1/3 0.2;0.25 2/3 0.4; ...
    0.75 1/9 0.6;0.125 4/9 0.8];
verifyEqual(testCase,points,expected,"AbsTol",1e-15);
end

function testQMCSequenceResumesFromPersistedSamplerState(testCase)
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()deleteStudyArtifacts(path));
study=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.QMCSampler(QMCType="sobol",Seed=9));
addNumericTrial(study);
first=collectNumericPoints(study,1);
verifyEqual(testCase,first,[0 0 0],"AbsTol",0);

resumed=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.QMCSampler(QMCType="sobol",Seed=9));
second=collectNumericPoints(resumed,1);
verifyEqual(testCase,second,[0.5 0.5 0.5],"AbsTol",0);
clear cleanup
end

function testQueuedTrialSurvivesStorageReload(testCase)
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()deleteStudyArtifacts(path));
study=radia.optuna.Study(StoragePath=path);
study.enqueueTrial(struct("x",0.375), ...
    UserAttrs=struct("source","baseline"));
reloaded=radia.optuna.Study(StoragePath=path);
trial=reloaded.ask();
verifyEqual(testCase,trial.suggestFloat("x",0,1),0.375);
verifyEqual(testCase,string(trial.UserAttrs.source),"baseline");
clear cleanup
end

function value=gridObjective(trial)
x=trial.suggest_int("x",1,2);
mode=trial.suggestCategorical("mode",{"a","b"});
value=x+double(string(mode)=="b");
end

function addNumericTrial(study)
trial=study.ask();
values=[trial.suggestFloat("x",0,1), ...
    trial.suggestFloat("y",0,1),trial.suggestFloat("z",0,1)];
study.tell(trial,sum(values));
end

function points=collectNumericPoints(study,count)
points=zeros(count,3);
for index=1:count
    trial=study.ask();
    points(index,:)=[trial.suggestFloat("x",0,1), ...
        trial.suggestFloat("y",0,1),trial.suggestFloat("z",0,1)];
    study.tell(trial,sum(points(index,:)));
end
end

function deleteStudyArtifacts(path)
for candidate=[path,path+".bak"]
    if isfile(candidate)
        delete(candidate);
    end
end
end
