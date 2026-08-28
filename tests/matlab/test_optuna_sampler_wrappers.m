function tests=test_optuna_sampler_wrappers
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
testCase.TestData.RemovePath= ...
    ~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemovePath, addpath(matlabDirectory); end
testCase.TestData.MatlabDirectory=matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testGPProposalResumesWithHyperparametersAndRandomState(testCase)
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()deleteStudyArtifacts(path));
continuous=radia.optuna.Study(Sampler=radia.optuna.GPSampler( ...
    Seed=31,NStartupTrials=3,CandidateCount=24,LocalSearchCount=1, ...
    DeterministicObjective=true,Backend="matlab-native"),AutoSave=false);
persisted=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.GPSampler(Seed=31,NStartupTrials=3,CandidateCount=24, ...
    LocalSearchCount=1,DeterministicObjective=true, ...
    Backend="matlab-native"));
for row=1:5
    first=continuous.ask();
    x1=first.suggestFloat("x",-1,1);
    continuous.tell(first,(x1-0.2)^2);
    second=persisted.ask();
    x2=second.suggestFloat("x",-1,1);
    persisted.tell(second,(x2-0.2)^2);
    verifyEqual(testCase,x1,x2,"AbsTol",0);
end
persisted=[]; %#ok<NASGU>
resumed=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.GPSampler(Seed=31,NStartupTrials=3,CandidateCount=24, ...
    LocalSearchCount=1,DeterministicObjective=true, ...
    Backend="matlab-native"));
first=continuous.ask();
x1=first.suggestFloat("x",-1,1);
second=resumed.ask();
x2=second.suggestFloat("x",-1,1);
verifyEqual(testCase,x1,x2,"AbsTol",0);
clear cleanup
end

function testAutoSamplerPolicyUsesGPAndNSGAIIIAtOptunaBudgets(testCase)
spec=struct("fixed_numeric",true,"dimensions",3, ...
    "has_constraints",false,"constraints_declared",true, ...
    "has_categorical",false,"is_conditional",false);
[name,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,1,100);
verifyEqual(testCase,name,"gp");
[name,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,1,500);
verifyEqual(testCase,name,"cmaes");
[name,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,3,100);
verifyEqual(testCase,name,"gp");
[name,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,3,500);
verifyEqual(testCase,name,"nsgaii");
[name,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,4,500);
verifyEqual(testCase,name,"nsgaiii");
spec.has_categorical=true;
spec.fixed_numeric=false;
[single,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,1,100);
[multi,~]=radia.optuna.internal.AutoSamplerPolicy.choose(spec,3,100);
verifyEqual(testCase,single,"tpe");
verifyEqual(testCase,multi,"motpe");
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
