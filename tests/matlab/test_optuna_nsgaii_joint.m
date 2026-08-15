function tests = test_optuna_nsgaii_joint
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(repositoryRoot,"matlab");
entries = string(strsplit(path,pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testConstraintDominationOrder(testCase)
study = radia.optuna.Study( ...
    Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler( ...
    Seed=41,PopulationSize=20),AutoSave=false);
values = [1 3;3 1;4 4;0 0;0.5 0.5;2 2;2.5 2.5];
constraints = {-Inf,0,-0.5,0.1,0.2,[],Inf};
for index = 1:size(values,1)
    trial = study.ask();
    study.tell(trial,values(index,:));
    if ~isempty(constraints{index})
        study.recordConstraints(trial,constraints{index});
    end
end

[feasible,violation,rank,crowding,order,missing] = ...
    radia.optuna.internal.ParetoSupport. ...
    constrainedRankAndCrowding(study,(0:6)',values);

verifyEqual(testCase,feasible, ...
    [true;true;true;false;false;false;false]);
verifyEqual(testCase,violation,[0;0;0;0.1;0.2;Inf;Inf]);
verifyEqual(testCase,rank(1:3),[1;1;2]);
verifyEqual(testCase,rank([4,5,7]),[3;4;5]);
verifyEqual(testCase,rank(6),6);
verifyTrue(testCase,all(isinf(crowding(1:3))));
verifyEqual(testCase,missing, ...
    [false;false;false;false;false;true;false]);
verifyEqual(testCase,order,[1;2;3;4;5;7;6]);
end

function testGenerationStateAndParentCacheAreFixed(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=47,PopulationSize=2, ...
    CrossoverProbability=1,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
for index = 1:2
    trial = study.ask();
    trial.suggestFloat("x",0,1);
    trial.suggestInteger("mesh",1,5);
    study.tell(trial,[index,3-index]);
end

firstChild = study.ask();
firstChild.suggestFloat("x",0,1);
firstChild.suggestInteger("mesh",1,5);
secondChild = study.ask();
secondChild.suggestFloat("x",0,1);
secondChild.suggestInteger("mesh",1,5);
state = study.SamplerStateTable.State{1};

verifyEqual(testCase,study.SamplerStateTable.Schema, ...
    "radia.optuna.nsgaii-sampler-state.v3");
verifyEqual(testCase,state.generation_by_trial, ...
    [0 0;1 0;2 1;3 1]);
verifyEqual(testCase,numel(state.generation_parent_cache),1);
verifyEqual(testCase,state.generation_parent_cache.generation,1);
verifyEqual(testCase,sort( ...
    state.generation_parent_cache.trial_numbers),[0;1]);
verifyTrue(testCase,all(ismember( ...
    firstChild.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));
verifyTrue(testCase,all(ismember( ...
    secondChild.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));

study.tell(firstChild,[0.5,0.5]);
thirdChild = study.ask();
thirdChild.suggestFloat("x",0,1);
thirdChild.suggestInteger("mesh",1,5);
verifyEqual(testCase,thirdChild.SystemAttrs.nsgaii_generation,1);
verifyTrue(testCase,all(ismember( ...
    thirdChild.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));
study.fail(secondChild,"test cleanup");
study.fail(thirdChild,"test cleanup");
end

function testParallelOvershootAndLateCompletionKeepParentCachesFixed(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=50,PopulationSize=2, ...
    CrossoverProbability=0,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);

initial = cell(4,1);
for index = 1:4
    initial{index} = study.ask();
    initial{index}.suggestFloat("x",0,1);
    verifyEqual(testCase, ...
        initial{index}.SystemAttrs.nsgaii_generation,0);
end
study.tell(initial{1},[0,3]);
study.tell(initial{2},[3,0]);

first = study.ask();
first.suggestFloat("x",0,1);
second = study.ask();
second.suggestFloat("x",0,1);
slow = study.ask();
slow.suggestFloat("x",0,1);
verifyEqual(testCase,[first.SystemAttrs.nsgaii_generation, ...
    second.SystemAttrs.nsgaii_generation, ...
    slow.SystemAttrs.nsgaii_generation],[1,1,1]);
verifyTrue(testCase,all(ismember( ...
    first.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));

% These late generation-zero results are deliberately better, but the
% already-created generation-one parent cache must remain unchanged.
study.tell(initial{3},[-100,-100]);
study.tell(initial{4},[-99,-99]);
lateCheck = study.ask();
lateCheck.suggestFloat("x",0,1);
verifyEqual(testCase,lateCheck.SystemAttrs.nsgaii_generation,1);
verifyTrue(testCase,all(ismember( ...
    lateCheck.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));
study.fail(lateCheck,"test cleanup");

study.tell(first,[1,2]);
study.tell(second,[2,1]);
generationTwo = study.ask();
generationTwo.suggestFloat("x",0,1);
verifyEqual(testCase,generationTwo.SystemAttrs.nsgaii_generation,2);
state = study.SamplerStateTable.State{1};
cacheIndex = find([state.generation_parent_cache.generation] == 2,1);
parentsBefore = sort( ...
    state.generation_parent_cache(cacheIndex).trial_numbers);
verifyFalse(testCase,ismember(slow.Number,parentsBefore));

study.tell(slow,[-200,-200]);
sameGeneration = study.ask();
sameGeneration.suggestFloat("x",0,1);
state = study.SamplerStateTable.State{1};
cacheIndex = find([state.generation_parent_cache.generation] == 2,1);
parentsAfter = sort( ...
    state.generation_parent_cache(cacheIndex).trial_numbers);
verifyEqual(testCase,parentsAfter,parentsBefore);
verifyFalse(testCase,ismember(slow.Number,parentsAfter));
study.fail(generationTwo,"test cleanup");
study.fail(sameGeneration,"test cleanup");
end

function testGenerationCacheAndRngResumeDeterministically(testCase)
primaryPath = string(tempname("C:\temp")) + ".mat";
clonePath = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@()cleanupStorage([primaryPath,clonePath]));

sampler = radia.optuna.NSGAIISampler(Seed=51,PopulationSize=2, ...
    CrossoverProbability=1,MutationProbability=0);
continuous = radia.optuna.Study( ...
    Name="nsgaii-resume",Directions=["minimize","minimize"], ...
    Sampler=sampler,StoragePath=primaryPath,AutoSave=true);
for index = 1:2
    trial = continuous.ask();
    trial.suggestFloat("x",-2,2);
    trial.suggestCategorical("material",["air","steel","magnet"]);
    continuous.tell(trial,[index,3-index]);
end
anchor = continuous.ask();
anchor.suggestFloat("x",-2,2);
anchor.suggestCategorical("material",["air","steel","magnet"]);
copyfile(primaryPath,clonePath,"f");

expected = continuous.ask();
expectedX = expected.suggestFloat("x",-2,2);
expectedMaterial = expected.suggestCategorical( ...
    "material",["air","steel","magnet"]);

resumed = radia.optuna.Study( ...
    Name="nsgaii-resume",Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler( ...
    Seed=51,PopulationSize=2,CrossoverProbability=1, ...
    MutationProbability=0),StoragePath=clonePath,AutoSave=true);
actual = resumed.ask();
actualX = actual.suggestFloat("x",-2,2);
actualMaterial = actual.suggestCategorical( ...
    "material",["air","steel","magnet"]);

verifyEqual(testCase,actual.Number,expected.Number);
verifyEqual(testCase,actual.SystemAttrs.nsgaii_generation, ...
    expected.SystemAttrs.nsgaii_generation);
verifyEqual(testCase,actual.SystemAttrs.nsgaii_parent_trial_numbers, ...
    expected.SystemAttrs.nsgaii_parent_trial_numbers);
verifyEqual(testCase,actualX,expectedX,AbsTol=0);
verifyEqual(testCase,string(actualMaterial),string(expectedMaterial));
clear cleanup
cleanupStorage([primaryPath,clonePath]);
end

function testMutationOmitsRelativeParametersForRandomFallback(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=48,PopulationSize=2, ...
    CrossoverProbability=1,MutationProbability=1);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
for index = 1:2
    trial = study.ask();
    trial.suggestFloat("x",-1,1);
    trial.suggestCategorical("material",["a","b"]);
    study.tell(trial,[index,3-index]);
end

child = study.ask();
verifyEqual(testCase,sort(child.SystemAttrs.nsgaii_mutated_parameters), ...
    ["material","x"]);
verifyFalse(testCase,isfield( ...
    child.SystemAttrs,"nsgaii_relative_search_space"));
x = child.suggestFloat("x",-1,1);
material = child.suggestCategorical("material",["a","b"]);
verifyGreaterThanOrEqual(testCase,x,-1);
verifyLessThanOrEqual(testCase,x,1);
verifyTrue(testCase,ismember(material,["a","b"]));
end

function testConstraintVectorShapeIsStable(testCase)
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler( ...
    Seed=49,PopulationSize=2),AutoSave=false);
first = study.ask();
study.tell(first,[0,1]);
study.recordConstraints(first,[0,0]);
second = study.ask();
study.tell(second,[1,0]);
study.recordConstraints(second,0);
verifyError(testCase,@() ...
    radia.optuna.internal.ParetoSupport.constrainedRankAndCrowding( ...
    study,[first.Number;second.Number],[0,1;1,0]), ...
    "radia:optuna:ConstraintShape");
end

function testJointMixedOffspringClonesOneWholeParent(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=42,PopulationSize=2, ...
    CrossoverProbability=0,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
parents = cell(2,3);
for index = 1:2
    trial = study.ask();
    parents{index,1} = trial.suggestFloat("x",-2,2);
    parents{index,2} = trial.suggest_int("mesh",1,7,Step=2);
    parents{index,3} = trial.suggestCategorical( ...
        "material",["air","steel","magnet"]);
    study.tell(trial,[index,3-index]);
end

child = study.ask();
offspring = {child.suggestFloat("x",-2,2), ...
    child.suggest_int("mesh",1,7,Step=2), ...
    child.suggestCategorical("material",["air","steel","magnet"])};
matches = false(2,1);
for index = 1:2
    matches(index) = offspring{1} == parents{index,1} && ...
        offspring{2} == parents{index,2} && ...
        string(offspring{3}) == string(parents{index,3});
end

verifyTrue(testCase,any(matches));
verifyEqual(testCase,child.SystemAttrs.nsgaii_sampling_mode,"joint");
verifyEqual(testCase,sort(child.SystemAttrs.nsgaii_relative_search_space), ...
    ["material","mesh","x"]);
verifyTrue(testCase,all(ismember( ...
    child.SystemAttrs.nsgaii_parent_trial_numbers,[0,1])));
verifyEqual(testCase,numel( ...
    child.SystemAttrs.nsgaii_parent_trial_numbers),1);
end

function testFeasibleParentBeatsInfeasibleObjectiveWinner(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=43,PopulationSize=2, ...
    CrossoverProbability=1,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
feasibleTrial = study.ask();
infeasibleTrial = study.ask();
moreInfeasibleTrial = study.ask();
feasibleX = feasibleTrial.suggestFloat("x",0,1);
infeasibleTrial.suggestFloat("x",0,1);
moreInfeasibleTrial.suggestFloat("x",0,1);
study.tell(feasibleTrial,[10,10]);
study.recordConstraints(feasibleTrial,0);
study.tell(infeasibleTrial,[0,0]);
study.recordConstraints(infeasibleTrial,0.1);
study.tell(moreInfeasibleTrial,[-1,-1]);
study.recordConstraints(moreInfeasibleTrial,0.2);

child = study.ask();
childX = child.suggestFloat("x",0,1); %#ok<NASGU>
verifyTrue(testCase,ismember(feasibleTrial.Number, ...
    child.SystemAttrs.nsgaii_parent_trial_numbers));
verifyFalse(testCase,ismember(moreInfeasibleTrial.Number, ...
    child.SystemAttrs.nsgaii_parent_trial_numbers));
verifyTrue(testCase,isfinite(feasibleX));
end

function testEmptyConstraintVectorIsPresentAndFeasible(testCase)
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler(Seed=91,PopulationSize=2), ...
    AutoSave=false);
emptyTrial = study.ask();
missingTrial = study.ask();
study.tell(emptyTrial,[1,1]);
study.recordConstraints(emptyTrial,zeros(1,0));
study.tell(missingTrial,[0,0]);
[present,values] = study.constraintRecord(emptyTrial.Number);
[missingPresent,~] = study.constraintRecord(missingTrial.Number);
[feasible,~,~,~,order,missing] = ...
    radia.optuna.internal.ParetoSupport.constrainedRankAndCrowding( ...
    study,[emptyTrial.Number;missingTrial.Number],[1,1;0,0]);
verifyTrue(testCase,present);
verifyEmpty(testCase,values);
verifyFalse(testCase,missingPresent);
verifyEqual(testCase,feasible,[true;false]);
verifyEqual(testCase,missing,[false;true]);
verifyEqual(testCase,order,[1;2]);
end

function testDynamicParameterFallsBackOutsideJointIntersection(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=44,PopulationSize=2, ...
    CrossoverProbability=0,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);

first = study.ask();
firstX = first.suggestFloat("x",-1,1);
first.suggestFloat("conditional",0,2);
study.tell(first,[0,1]);
second = study.ask();
secondX = second.suggestFloat("x",-1,1);
study.tell(second,[1,0]);

child = study.ask();
childX = child.suggestFloat("x",-1,1);
conditional = child.suggestFloat("conditional",0,2);
verifyEqual(testCase,child.SystemAttrs.nsgaii_relative_search_space,"x");
verifyTrue(testCase,childX == firstX || childX == secondX);
verifyGreaterThanOrEqual(testCase,conditional,0);
verifyLessThanOrEqual(testCase,conditional,2);
end

function testRunningAndFailedTrialsAreNeverParents(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=45,PopulationSize=2, ...
    CrossoverProbability=0,MutationProbability=0);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
for index = 1:2
    trial = study.ask();
    trial.suggestFloat("x",0,1);
    study.tell(trial,[index,3-index]);
end
running = study.ask();
running.suggestFloat("x",0,1);
failed = study.ask();
failed.suggestFloat("x",0,1);
study.fail(failed,"worker failed");

child = study.ask();
child.suggestFloat("x",0,1);
parents = child.SystemAttrs.nsgaii_parent_trial_numbers;
verifyTrue(testCase,all(ismember(parents,[0,1])));
verifyFalse(testCase,any(ismember(parents,[running.Number,failed.Number])));
study.fail(running,"test cleanup");
end

function testConstraintCallbackPersistsPerCompletedTrial(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=46,PopulationSize=2, ...
    ConstraintsFcn=@(trial)trial.Params.x-0.4);
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
study.optimize(@constrainedObjective,4);

verifyEqual(testCase,height(study.ConstraintTable),4);
verifyTrue(testCase,all(isfinite(study.ConstraintTable.Value)));
verifyEqual(testCase,study.ConstraintTable.TrialNumber,(0:3)');
end

function testOptuna49DefaultsAndCrossoverSurface(testCase)
sampler = radia.optuna.NSGAIISampler();
verifyEqual(testCase,sampler.PopulationSize,50);
verifyTrue(testCase,isnan(sampler.MutationProbability));
verifyEqual(testCase,sampler.CrossoverProbability,0.9);
verifyEqual(testCase,sampler.SwappingProbability,0.5);
verifyClass(testCase,sampler.Crossover, ...
    "radia.optuna.nsgaii.UniformCrossover");
crossovers = {radia.optuna.nsgaii.UniformCrossover(), ...
    radia.optuna.nsgaii.BLXAlphaCrossover(), ...
    radia.optuna.nsgaii.SPXCrossover(), ...
    radia.optuna.nsgaii.SBXCrossover(), ...
    radia.optuna.nsgaii.VSBXCrossover(), ...
    radia.optuna.nsgaii.UNDXCrossover()};
verifyEqual(testCase,cellfun(@(item)item.NParents,crossovers), ...
    [2,2,3,2,2,3]);
verifyError(testCase,@()radia.optuna.NSGAIISampler( ...
    PopulationSize=2,Crossover=radia.optuna.nsgaii.SPXCrossover()), ...
    "radia:optuna:NSGAIIPopulation");
end

function testAllCrossoversProduceContainedMixedChildren(testCase)
crossovers = {radia.optuna.nsgaii.UniformCrossover(), ...
    radia.optuna.nsgaii.BLXAlphaCrossover(), ...
    radia.optuna.nsgaii.SPXCrossover(), ...
    radia.optuna.nsgaii.SBXCrossover(), ...
    radia.optuna.nsgaii.VSBXCrossover(), ...
    radia.optuna.nsgaii.UNDXCrossover()};
for crossoverIndex = 1:numel(crossovers)
    crossover = crossovers{crossoverIndex};
    populationSize = max(3,crossover.NParents);
    sampler = radia.optuna.NSGAIISampler(Seed=100+crossoverIndex, ...
        PopulationSize=populationSize,Crossover=crossover, ...
        CrossoverProbability=1,MutationProbability=0);
    study = radia.optuna.Study(Directions=["minimize","minimize"], ...
        Sampler=sampler,AutoSave=false);
    for index = 1:populationSize
        trial = study.ask();
        x = trial.suggestFloat("x",1e-3,1e3,Log=true);
        mesh = trial.suggest_int("mesh",1,9,Step=2);
        mode = trial.suggestCategorical("mode",["A","B"]);
        study.tell(trial,[x,(mesh-5)^2+double(mode=="B")]);
    end
    child = study.ask();
    x = child.suggestFloat("x",1e-3,1e3,Log=true);
    mesh = child.suggest_int("mesh",1,9,Step=2);
    mode = child.suggestCategorical("mode",["A","B"]);
    verifyGreaterThanOrEqual(testCase,x,1e-3);
    verifyLessThanOrEqual(testCase,x,1e3);
    verifyTrue(testCase,ismember(mesh,1:2:9));
    verifyTrue(testCase,ismember(mode,["A","B"]));
    verifyEqual(testCase,numel( ...
        child.SystemAttrs.nsgaii_parent_trial_numbers), ...
        crossover.NParents);
end
end

function testPrunedTrialRunsConstraintCallback(testCase)
sampler = radia.optuna.NSGAIISampler(Seed=121,PopulationSize=2, ...
    ConstraintsFcn=@(~)zeros(1,0));
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
trial = study.ask();
trial.prune();
[present,values] = study.constraintRecord(trial.Number);
verifyTrue(testCase,present);
verifyEmpty(testCase,values);
verifyEqual(testCase,study.TrialTable.State,"PRUNED");
end

function testEliteSelectionUsesCrowdingOnlyAtCutoff(testCase)
study = radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler(Seed=122,PopulationSize=3), ...
    AutoSave=false);
values = [0,4;1,3;2,2;3,1;4,0];
trials = zeros(5,1);
for index = 1:5
    trial = study.ask();
    trials(index) = trial.Number;
    study.tell(trial,values(index,:));
end
order = radia.optuna.internal.ParetoSupport.eliteSelectionOrder( ...
    study,trials,values,3);
verifyEqual(testCase,numel(order),3);
verifyEqual(testCase,order,[1;5;2]);
verifyEqual(testCase,numel(unique(order)),3);
end

function values = constrainedObjective(trial)
x = trial.suggestFloat("x",0,1);
values = [x^2,(x-1)^2];
end

function cleanupStorage(paths)
for path = reshape(string(paths),1,[])
    for candidate = [path,path + ".bak"]
        if isfile(candidate)
            delete(candidate);
        end
    end
end
end
