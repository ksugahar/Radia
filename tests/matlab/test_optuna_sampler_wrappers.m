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

function testBruteForceExhaustsConditionalDefineByRunTree(testCase)
sampler=radia.optuna.BruteForceSampler(Seed=11);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
result=study.optimize(@conditionalFiniteObjective,100);
verifyEqual(testCase,height(result),6);
verifyTrue(testCase,all(result.State=="COMPLETE"));
keys=strings(1,height(result));
for row=1:height(result)
    params=result.Params{row};
    if string(params.kind)=="float"
        keys(row)="float:"+string(params.x);
    else
        keys(row)="int:"+string(params.a)+":"+string(params.b);
    end
end
verifyEqual(testCase,numel(unique(keys)),6);
expected=sort(["float:1","float:1.5","float:2", ...
    "int:1:1","int:1:2","int:2:2"]);
verifyEqual(testCase,sort(keys),expected);
end

function testBruteForceAvoidsRunningDuplicateAndRejectsInfiniteSpace(testCase)
study=radia.optuna.Study(Sampler= ...
    radia.optuna.BruteForceSampler(Seed=5),AutoSave=false);
first=study.ask();
a=first.suggestCategorical("choice",[1 2]);
second=study.ask();
b=second.suggestCategorical("choice",[1 2]);
verifyNotEqual(testCase,a,b);
third=study.ask();
verifyError(testCase,@()third.suggestFloat("x",0,1), ...
    "radia:optuna:BruteForceInfinite");
end

function testNSGAIIIReferenceNichingAndGenerationMetadata(testCase)
references=[1 0 0;0 1 0;0 0 1;1 1 1];
sampler=radia.optuna.NSGAIIISampler(Seed=17,PopulationSize=4, ...
    ReferencePoints=references,MutationProbability=0);
study=radia.optuna.Study(Sampler=sampler, ...
    Directions=["minimize","minimize","minimize"],AutoSave=false);
values=[10 0 0;0 10 0;0 0 10;3 3 3;9 0.1 0.1;0.1 9 0.1];
numbers=zeros(size(values,1),1);
for row=1:size(values,1)
    violation=-1;
    if row==1, violation=1; end
    frozen=radia.optuna.createTrial(Values=values(row,:),State="COMPLETE", ...
        Constraints=violation,ConstraintPresent=true);
    addTrial(study,frozen);
    numbers(row)=study.TrialTable.TrialNumber(end);
end

selected=sampler.selectElitePopulation(study,numbers);
verifyEqual(testCase,numel(selected),4);
verifyEqual(testCase,numel(unique(selected)),4);
verifyTrue(testCase,ismember(numbers(4),selected));
verifyFalse(testCase,ismember(numbers(1),selected));

liveStudy=radia.optuna.Study(Sampler=radia.optuna.NSGAIIISampler( ...
    Seed=3,PopulationSize=2),Directions=["minimize","minimize"], ...
    AutoSave=false);
for row=1:2
    trial=liveStudy.ask();
    trial.suggestFloat("x",0,1);
    liveStudy.tell(trial,[row,3-row]);
end
child=liveStudy.ask();
child.suggestFloat("x",0,1);
verifyEqual(testCase,child.SystemAttrs.nsgaiii_generation,1);
verifyEqual(testCase,child.SystemAttrs.nsgaiii_elite_strategy, ...
    "reference_line_niching");
end

function testNSGAIIIWorseFrontCannotMoveCutoffNormalization(testCase)
references=[1 0 0;0 1 0;0 0 1;1 1 1];
base=[0 8 8;8 0 8;8 8 0;3 3 3;2 6 6;6 2 6];
first=nsga3Selection(base,references,29);
second=nsga3Selection([base;1e12 1e12 1e12],references,29);
verifyEqual(testCase,second,first);
end

function testNSGAIIITransformsMatchOptuna49Fixtures(testCase)
references=radia.optuna.internal.NSGAIIISupport. ...
    defaultReferencePoints(3,3);
verifyEqual(testCase,references,[3 0 0;2 1 0;2 0 1;1 2 0; ...
    1 1 1;1 0 2;0 3 0;0 2 1;0 1 2;0 0 3]);

values=[0 0 4;0 4 0;4 0 0;1 1 1;2 2 2];
normalized=radia.optuna.internal.NSGAIIISupport. ...
    normalizeObjectives(values);
verifyEqual(testCase,normalized,[0 0 1;0 1 0;1 0 0; ...
    .25 .25 .25;.5 .5 .5],"AbsTol",1e-15);
[associations,distances]=radia.optuna.internal.NSGAIIISupport. ...
    associate(normalized,[1 0 0;0 1 0;0 0 1;1 1 1]);
verifyEqual(testCase,associations,[3;2;1;4;4]);
verifyEqual(testCase,distances,zeros(5,1),"AbsTol",1e-15);

nonfinite=[-Inf 5 5;5 -Inf 5;5 5 -Inf;2 2 2;Inf 1 1;1 Inf 1];
actual=radia.optuna.internal.NSGAIIISupport. ...
    normalizeObjectives(nonfinite);
expected=[0 .07692307692307682 1.076923076923077; ...
    .0769230769230769 0 1.076923076923077; ...
    .0769230769230769 .07692307692307682 0; ...
    .0625 .0625 .875; ...
    .13461538461538458 .05769230769230761 .8076923076923078; ...
    .057692307692307675 .13461538461538441 .8076923076923078];
verifyEqual(testCase,actual,expected,"AbsTol",2e-15);
end

function testGPMaternConstraintsPendingAndMultiobjective(testCase)
sampler=radia.optuna.GPSampler(Seed=23,NStartupTrials=3, ...
    CandidateCount=24,LocalSearchCount=1,DeterministicObjective=true, ...
    ConstraintsFcn=@(trial)0.2-trial.Params.x);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
for row=1:4
    trial=study.ask();
    x=trial.suggestFloat("x",0,1);
    study.tell(trial,(x-0.1)^2);
end
first=study.ask();
x1=first.suggestFloat("x",0,1);
second=study.ask();
x2=second.suggestFloat("x",0,1);
verifyNotEqual(testCase,x1,x2);
verifyEqual(testCase,height(study.ConstraintTable),4);
verifyTrue(testCase,any(study.SystemAttrTable.Name=="gp_acquisition"));

multi=radia.optuna.Study(Sampler=radia.optuna.GPSampler( ...
    Seed=7,NStartupTrials=3,CandidateCount=20,LocalSearchCount=1, ...
    MonteCarloSamples=4,DeterministicObjective=true), ...
    Directions=["minimize","minimize"],AutoSave=false);
for row=1:4
    trial=multi.ask();
    x=trial.suggestFloat("x",0,1);
    multi.tell(trial,[x^2,(1-x)^2]);
end
rows=multi.SystemAttrTable.Name=="gp_acquisition";
verifyTrue(testCase,any(rows));
verifyEqual(testCase,string(jsondecode( ...
    multi.SystemAttrTable.ValueJSON(find(rows,1,"last")))), ...
    "expected_hypervolume_improvement");
end

function testGPCholeskyJitterKeepsDuplicateObservationsFinite(testCase)
x=[0;0;0.5;1];
y=[1;1;0.2;0.8];
[model,theta]=radia.optuna.internal.GaussianProcess.fit( ...
    x,y,false,true);
[meanValue,stdValue]=radia.optuna.internal.GaussianProcess.predict( ...
    model,[0;0.25;1]);
verifyTrue(testCase,all(isfinite(meanValue)));
verifyTrue(testCase,all(isfinite(stdValue)));
verifyGreaterThan(testCase,min(stdValue),0);
verifyTrue(testCase,all(isfinite(theta)));
end

function testGPMatchesOptuna49KernelHyperparameters(testCase)
qmc=radia.optuna.QMCSampler(QMCType="sobol");
x=qmc.unitPoints(2,10);
y=sin(5*x(:,1))+0.3*(x(:,2)-0.4).^2;
[~,theta]=radia.optuna.internal.GaussianProcess.fit( ...
    x,y,[false false],true);
lengthScales=reshape(1./sqrt(exp(theta(1:2))),1,[]);
kernelScale=exp(theta(3));
% Optuna 4.9 official GP fit on the identical standardized fixture.
verifyEqual(testCase,lengthScales,[0.450364137125575,3.69198874495208], ...
    "RelTol",5e-3);
verifyEqual(testCase,kernelScale,1.32415059683817,"RelTol",5e-3);
end

function testGPHartmann6FixedBudgetQualityMatchesOptuna49(testCase)
startup=10;
sampler=radia.optuna.GPSampler(Seed=1,NStartupTrials=startup, ...
    DeterministicObjective=true);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
points=radia.optuna.QMCSampler(QMCType="sobol").unitPoints(6,startup);
for row=1:25
    if row<=startup
        params=struct;
        for dimension=1:6
            params.("x"+dimension)=points(row,dimension);
        end
        study.enqueueTrial(params);
    end
    trial=study.ask();
    values=zeros(1,6);
    for dimension=1:6
        values(dimension)=trial.suggestFloat("x"+dimension,0,1);
    end
    study.tell(trial,hartmann6(values));
end
regret=study.bestValue()-(-3.322368011415515);
% The official Optuna 4.9 result is 0.4034 for this fixed-start budget;
% allow margin for the native acquisition optimizer while rejecting the
% pre-parity regression (>1.5).
verifyLessThan(testCase,regret,0.8);
end

function testGPThreeObjectiveFixedBudgetQualityMatchesOptuna49(testCase)
startup=10;
reference=[2.1 2.1 2.1];
sampler=radia.optuna.GPSampler(Seed=1,NStartupTrials=startup, ...
    DeterministicObjective=true);
study=radia.optuna.Study(Sampler=sampler, ...
    Directions=["minimize","minimize","minimize"],AutoSave=false);
points=radia.optuna.QMCSampler(QMCType="sobol").unitPoints(2,startup);
for row=1:25
    if row<=startup
        study.enqueueTrial(struct("x",points(row,1),"y",points(row,2)));
    end
    trial=study.ask();
    x=trial.suggestFloat("x",0,1);
    y=trial.suggestFloat("y",0,1);
    study.tell(trial,[x^2+y^2,(x-1)^2+y^2,x^2+(y-1)^2]);
end
front=study.paretoFront();
values=zeros(height(front),3);
for row=1:height(front)
    for objective=1:3
        mask=study.ObjectiveTable.TrialNumber==front.TrialNumber(row) & ...
            study.ObjectiveTable.ObjectiveIndex==objective;
        values(row,objective)=study.ObjectiveTable.Value(mask);
    end
end
hypervolume=radia.optuna.internal.ParetoSupport. ...
    computeHypervolume(values,reference);
% Optuna 4.9 reaches 6.9198 and the native sampler 6.8772 on this exact
% 10-point Sobol startup plus 15 adaptive-trial fixture.
verifyGreaterThan(testCase,hypervolume,6.8);
end

function testFastTwoObjectiveHypervolumeImprovementIsExact(testCase)
front=[0.2 0.8;0.5 0.4;0.8 0.2];
reference=[1 1];
samples=[0.1 0.9;0.3 0.3;0.9 0.1;1.1 0.1;0.6 0.7];
actual=radia.optuna.internal.ParetoSupport. ...
    hypervolumeImprovement2D(samples,front,reference);
base=radia.optuna.internal.ParetoSupport. ...
    computeHypervolume(front,reference);
expected=zeros(size(samples,1),1);
for row=1:size(samples,1)
    expected(row)=max(0,radia.optuna.internal.ParetoSupport. ...
        computeHypervolume([front;samples(row,:)],reference)-base);
end
verifyEqual(testCase,actual,expected,"AbsTol",1e-14);
end

function testGPProposalResumesWithHyperparametersAndRandomState(testCase)
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()deleteStudyArtifacts(path));
continuous=radia.optuna.Study(Sampler=radia.optuna.GPSampler( ...
    Seed=31,NStartupTrials=3,CandidateCount=24,LocalSearchCount=1, ...
    DeterministicObjective=true),AutoSave=false);
persisted=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.GPSampler(Seed=31,NStartupTrials=3,CandidateCount=24, ...
    LocalSearchCount=1,DeterministicObjective=true));
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
    LocalSearchCount=1,DeterministicObjective=true));
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

function testSobolSequenceMatchesOptunaUnscrambledPrefix(testCase)
sampler=radia.optuna.QMCSampler(QMCType="sobol",Seed=7);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
addNumericTrial(study); % Initial trial is sampled independently.
points=collectNumericPoints(study,5);
expected=[0 0 0;0.5 0.5 0.5;0.75 0.25 0.25; ...
    0.25 0.75 0.75;0.375 0.375 0.625];
verifyEqual(testCase,points,expected,"AbsTol",0);
end

function testSobolSequenceMatchesOptunaThroughDimension32(testCase)
sampler=radia.optuna.QMCSampler(QMCType="sobol");
points=sampler.unitPoints(32,64);
rows=[1 2 3 4 5 8 16 32 64];
columns=[1 2 3 4 8 16 24 32];
expected=[ ...
    0 0 0 0 0 0 0 0; ...
    .5 .5 .5 .5 .5 .5 .5 .5; ...
    .75 .25 .25 .25 .75 .25 .25 .25; ...
    .25 .75 .75 .75 .25 .75 .75 .75; ...
    .375 .375 .625 .875 .875 .875 .875 .125; ...
    .125 .625 .375 .125 .625 .125 .125 .875; ...
    .0625 .9375 .5625 .3125 .3125 .8125 .8125 .1875; ...
    .03125 .53125 .90625 .96875 .53125 .84375 .28125 .40625; ...
    .015625 .796875 .359375 .453125 .140625 .765625 .546875 .921875];
verifyEqual(testCase,points(rows,columns),expected,"AbsTol",0);
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

function testScrambledQMCIsSeededBoundedAndNondegenerate(testCase)
for qmcType=["sobol","halton"]
    first=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
        QMCType=qmcType,Scramble=true,Seed=123),AutoSave=false);
    second=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
        QMCType=qmcType,Scramble=true,Seed=123),AutoSave=false);
    different=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
        QMCType=qmcType,Scramble=true,Seed=124),AutoSave=false);
    addNumericTrial(first);
    addNumericTrial(second);
    addNumericTrial(different);
    a=collectNumericPoints(first,4);
    b=collectNumericPoints(second,4);
    c=collectNumericPoints(different,4);
    verifyEqual(testCase,a,b,"AbsTol",0);
    verifyNotEqual(testCase,a,c);
    verifyTrue(testCase,all(a>=0 & a<1,"all"));
    verifyGreaterThan(testCase,numel(unique(a(:,1))),1);
end
end

function testSobolScramblePreservesDigitalNetBalance(testCase)
study=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
    QMCType="sobol",Scramble=true,Seed=321),AutoSave=false);
trial=study.ask();
trial.suggestFloat("x",0,1);
trial.suggestFloat("y",0,1);
study.tell(trial,0);
counts=zeros(4);
for row=1:16
    trial=study.ask();
    x=trial.suggestFloat("x",0,1);
    y=trial.suggestFloat("y",0,1);
    study.tell(trial,0);
    ix=min(floor(4*x)+1,4);
    iy=min(floor(4*y)+1,4);
    counts(ix,iy)=counts(ix,iy)+1;
end
verifyEqual(testCase,counts,ones(4));
end

function testScrambledQMCDiscrepancyMatchesSciPyQualityBand(testCase)
sobol=radia.optuna.QMCSampler(QMCType="sobol",Scramble=true,Seed=321);
halton=radia.optuna.QMCSampler(QMCType="halton",Scramble=true,Seed=321);
sobolDiscrepancy=centeredDiscrepancy(sobol.unitPoints(5,256));
haltonDiscrepancy=centeredDiscrepancy(halton.unitPoints(5,256));
verifyGreaterThan(testCase,sobolDiscrepancy,2.0e-4);
verifyLessThan(testCase,sobolDiscrepancy,3.0e-4);
verifyGreaterThan(testCase,haltonDiscrepancy,3.0e-4);
verifyLessThan(testCase,haltonDiscrepancy,5.0e-4);
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

function value=conditionalFiniteObjective(trial)
kind=string(trial.suggestCategorical("kind",{"float","int"}));
if kind=="float"
    value=trial.suggestFloat("x",1,2,Step=0.5);
else
    a=trial.suggest_int("a",1,2);
    b=trial.suggest_int("b",a,2);
    value=a+b;
end
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

function selected=nsga3Selection(values,references,seed)
sampler=radia.optuna.NSGAIIISampler(Seed=seed,PopulationSize=4, ...
    ReferencePoints=references);
study=radia.optuna.Study(Sampler=sampler, ...
    Directions=["minimize","minimize","minimize"],AutoSave=false);
numbers=zeros(size(values,1),1);
for row=1:size(values,1)
    addTrial(study,radia.optuna.createTrial( ...
        Values=values(row,:),State="COMPLETE"));
    numbers(row)=study.TrialTable.TrialNumber(end);
end
selected=reshape(sampler.selectElitePopulation(study,numbers),1,[]);
end

function value=hartmann6(x)
alpha=[1 1.2 3 3.2];
a=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14; ...
    3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
p=1e-4*[1312 1696 5569 124 8283 5886; ...
    2329 4135 8307 3736 1004 9991; ...
    2348 1451 3522 2883 3047 6650; ...
    4047 8828 8732 5743 1091 381];
value=-sum(alpha.*exp(-sum(a.*(x-p).^2,2)'));
end

function value=centeredDiscrepancy(points)
[count,dimension]=size(points);
centered=abs(points-.5);
first=(13/12)^dimension;
second=(2/count)*sum(prod(1+.5*centered-.5*centered.^2,2));
pairSum=0;
for row=1:count
    factors=1+.5*centered(row,:)+.5*centered- ...
        .5*abs(points(row,:)-points);
    pairSum=pairSum+sum(prod(factors,2));
end
value=first-second+pairSum/count^2;
end

function deleteStudyArtifacts(path)
for candidate=[path,path+".bak"]
    if isfile(candidate)
        delete(candidate);
    end
end
end
