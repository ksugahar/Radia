function tests=test_optuna_upstream_oracle
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

function testOracleProvenance(testCase)
oracle=testCase.TestData.Oracle;
verifyEqual(testCase,string(oracle.schema), ...
    "radia.test.optuna-upstream-oracle.v1");
verifyEqual(testCase,string(oracle.optuna_version),"4.9.0");
verifyNotEmpty(testCase,string(oracle.numpy_version));
verifyNotEmpty(testCase,string(oracle.scipy_version));
verifyNotEmpty(testCase,string(oracle.python_version));
verifyNotEmpty(testCase,string(oracle.torch_version));
verifyNotEmpty(testCase,string(oracle.cmaes_version));
end

function testOraclePolicyManifestIsComplete(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
fixtureDirectory=fullfile(root,"tests","matlab","fixtures");
manifest=jsondecode(fileread(fullfile( ...
    fixtureDirectory,"optuna_test_manifest.json")));
verifyEqual(testCase,string(manifest.schema), ...
    "radia.test.optuna-matlab-policy.v1");
verifyEqual(testCase,string(manifest.upstream_version),"4.9.0");

files=dir(fullfile(root,"tests","matlab","test_optuna*.m"));
actual=strings(0,1);
for fileIndex=1:numel(files)
    source=fileread(fullfile(files(fileIndex).folder,files(fileIndex).name));
    names=regexp(source, ...
        '(?m)^function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(', ...
        'tokens');
    for nameIndex=1:numel(names)
        if string(names{nameIndex}{1})~=erase(string(files(fileIndex).name),".m")
            actual(end+1,1)=string(files(fileIndex).name)+"::"+ ...
                string(names{nameIndex}{1}); %#ok<AGROW>
        end
    end
end

entries=manifest.entries;
declared=strings(numel(entries),1);
allowed=["upstream-python","upstream-mcp","matlab-integration"];
for index=1:numel(entries)
    declared(index)=string(entries(index).file)+"::"+ ...
        string(entries(index).test);
    verifyTrue(testCase,ismember(string(entries(index).classification),allowed));
    if string(entries(index).classification)=="upstream-python"
        verifyEqual(testCase,string(entries(index).oracle), ...
            "optuna49_oracle.json");
    elseif string(entries(index).classification)=="upstream-mcp"
        verifyEqual(testCase,string(entries(index).oracle), ...
            "optuna49_mcp_oracle.json");
    else
        verifyTrue(testCase,isempty(entries(index).oracle));
        verifyNotEmpty(testCase,string(entries(index).scope));
    end
end
verifyEqual(testCase,sort(declared),sort(actual));
verifyEqual(testCase,numel(unique(declared)),numel(declared));
end

function testNumpyRandomStateSeedContract(testCase)
expected=testCase.TestData.Oracle.numpy_random_state_seed_contract;
stream=radia.optuna.internal.NumpyRandomState(37);
uniforms=rand(stream,400,1);
positions=double(expected.uniform_positions_zero_based)+1;
verifyEqual(testCase,uniforms(positions), ...
    reshape(double(expected.uniform_values),[],1),AbsTol=0);

stream=radia.optuna.internal.NumpyRandomState(37);
verifyEqual(testCase,randn(stream,7,1), ...
    reshape(double(expected.normal_values),[],1),AbsTol=0);
stream=radia.optuna.internal.NumpyRandomState(123);
verifyEqual(testCase,randperm(stream,10), ...
    reshape(double(expected.permutation_one_based),1,[]));
stream=radia.optuna.internal.NumpyRandomState(123);
verifyEqual(testCase,randi(stream,10,1,8), ...
    reshape(double(expected.integers_one_based),1,[]));
end

function testUpstreamDefaults(testCase)
expected=testCase.TestData.Oracle.defaults;
single=radia.optuna.create_study(AutoSave=false);
multi=radia.optuna.create_study( ...
    directions=["minimize","maximize"],AutoSave=false);
verifyTrue(testCase,startsWith(single.Name,string(expected.anonymous_prefix)));
verifyEqual(testCase,string(class(single.Sampler)), ...
    "radia.optuna."+string(expected.single_sampler));
verifyEqual(testCase,string(class(multi.Sampler)), ...
    "radia.optuna."+string(expected.multi_sampler));
verifyEqual(testCase,multi.Sampler.PopulationSize, ...
    double(expected.multi_population_size));
verifyEqual(testCase,string(class(single.Pruner)), ...
    "radia.optuna."+string(expected.pruner));
end

function testUnseededSamplerDefaultsMatchUpstreamSemantics(testCase)
expected=testCase.TestData.Oracle.sampler_seed_defaults;
names=["RandomSampler","TPESampler","CmaEsSampler","GPSampler", ...
    "GridSampler","NSGAIISampler","NSGAIIISampler","QMCSampler", ...
    "BruteForceSampler"];
for index=1:numel(names)
    entry=expected.constructors.(names(index));
    verifyTrue(testCase,logical(entry.default_is_none));
    verifyEqual(testCase,string(entry.parameter),"seed");
end
verifyTrue(testCase,logical(expected.exact_sequence_oracle_requires_explicit_seed));

globalStateBefore=rng;
first={radia.optuna.RandomSampler(),radia.optuna.TPESampler(), ...
    radia.optuna.CmaEsSampler(),radia.optuna.GPSampler(), ...
    radia.optuna.GridSampler(struct("x",{{0,1}})), ...
    radia.optuna.NSGAIISampler(),radia.optuna.NSGAIIISampler(), ...
    radia.optuna.QMCSampler(),radia.optuna.BruteForceSampler()};
second={radia.optuna.RandomSampler(),radia.optuna.TPESampler(), ...
    radia.optuna.CmaEsSampler(),radia.optuna.GPSampler(), ...
    radia.optuna.GridSampler(struct("x",{{0,1}})), ...
    radia.optuna.NSGAIISampler(),radia.optuna.NSGAIIISampler(), ...
    radia.optuna.QMCSampler(),radia.optuna.BruteForceSampler()};
globalStateAfter=rng;
firstSeeds=cellfun(@(sampler)sampler.Seed,first);
secondSeeds=cellfun(@(sampler)sampler.Seed,second);
verifyTrue(testCase,all(firstSeeds~=secondSeeds));
verifyEqual(testCase,globalStateAfter,globalStateBefore);

verifyEqual(testCase,radia.optuna.RandomSampler(37).Seed,37);
verifyEqual(testCase,radia.optuna.TPESampler(Seed=37).Seed,37);
end

function testUpstreamTellContract(testCase)
expected=testCase.TestData.Oracle.tell;
study=radia.optuna.Study(AutoSave=false);
trial=study.ask();
trial.report(3.5,0);
trial.report(2.5,1);
snapshot=study.tell(trial,State="PRUNED");
verifyEqual(testCase,snapshot.State,string(expected.pruned_state));
verifyEqual(testCase,snapshot.Value,double(expected.pruned_value));

trial=study.ask();
snapshot=study.tell(trial,State="FAIL");
verifyEqual(testCase,snapshot.State,string(expected.failed_state));
verifyEqual(testCase,isnan(snapshot.Value),logical(expected.failed_value_is_none));

trial=study.ask();
warning("off","radia:optuna:InvalidObjectiveValue");
cleanup=onCleanup(@()warning("on","radia:optuna:InvalidObjectiveValue"));
snapshot=study.tell(trial);
verifyEqual(testCase,snapshot.State,string(expected.missing_state));
clear cleanup

trial=study.ask();
snapshot=study.tell(trial,Inf);
verifyEqual(testCase,snapshot.State,string(expected.infinite_state));
verifyEqual(testCase,string(expected.infinite_value),"Infinity");
verifyEqual(testCase,snapshot.Value,Inf);

trial=study.ask();
trial.report(4.5,1);
snapshot=study.tell(trial.Number,4.25);
verifyEqual(testCase,snapshot.Number,double(expected.by_number.number));
verifyEqual(testCase,snapshot.State,string(expected.by_number.state));
verifyEqual(testCase,snapshot.Value,double(expected.by_number.value));
verifyEqual(testCase,snapshot.last_step(),double(expected.by_number.last_step));
end

function testTrialPrunedExceptionMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.trial_pruned_exception;
callbackRows=struct("state",{},"value",{},"last_step",{});
study=radia.optuna.Study(AutoSave=false);
study.optimize(@objective,1,Callbacks={@callback});
snapshot=study.trials("PRUNED");
verifyEqual(testCase,height(snapshot),1);
verifyEqual(testCase,snapshot.State,string(expected.state));
verifyEqual(testCase,snapshot.Value,double(expected.value));
frozen=study.get_trials("PRUNED");
verifyEqual(testCase,numel(frozen),1);
verifyEqual(testCase,frozen(1).last_step(),double(expected.last_step));
verifyEqual(testCase,numel(callbackRows),double(expected.callback_count));
verifyEqual(testCase,callbackRows(1).state,string(expected.callback.state));
verifyEqual(testCase,callbackRows(1).value,double(expected.callback.value));
verifyEqual(testCase,callbackRows(1).last_step, ...
    double(expected.callback.last_step));

    function value=objective(trial) %#ok<STOUT>
        trial.report(7,0);
        trial.report(3,2);
        throw(radia.optuna.TrialPruned("oracle prune"));
    end

    function callback(~,trial)
        callbackRows(end+1)=struct("state",trial.State, ...
            "value",trial.Value,"last_step",trial.last_step());
    end
end

function testFinishedTellAndConstraintFailureMatchUpstream(testCase)
expected=testCase.TestData.Oracle.lifecycle_errors;
study=radia.optuna.Study(AutoSave=false);
trial=study.ask();
study.tell(trial,1);
verifyError(testCase,@()study.tell(trial,2),"radia:optuna:TrialState");
snapshot=study.tell(trial,999,SkipIfFinished=true);
verifyEqual(testCase,snapshot.Value,double(expected.finished_tell.skip_value));
verifyEqual(testCase,snapshot.State,string(expected.finished_tell.skip_state));

sampler=radia.optuna.TPESampler(Seed=83, ...
    ConstraintsFcn=@(~)failingConstraintCallback());
constrained=radia.optuna.Study(Sampler=sampler,AutoSave=false);
verifyError(testCase,@()constrained.optimize( ...
    @(item)item.suggest_float("x",0,1),1), ...
    "radia:test:ConstraintCallback");
verifyEqual(testCase,constrained.TrialTable.State, ...
    string(expected.constraint_callback_failure.state));
verifyEqual(testCase,isfinite(constrained.TrialTable.Value), ...
    logical(expected.constraint_callback_failure.value_is_finite));
verifyEmpty(testCase,constrained.ConstraintTable);
end

function testFixedTrialMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.fixed_trial;
trial=radia.optuna.FixedTrial( ...
    struct("x",0.5,"n",3,"kind","b"),Number=7);
verifyEqual(testCase,trial.Number,double(expected.number));
verifyEqual(testCase,numel(fieldnames(trial.Params)), ...
    numel(fieldnames(expected.initial_params)));
verifyEqual(testCase,numel(fieldnames(trial.Distributions)), ...
    double(expected.initial_distribution_count));
verifyEqual(testCase,trial.suggest_float("x",0,1), ...
    double(expected.values.x));
verifyEqual(testCase,trial.suggest_int("n",1,5), ...
    double(expected.values.n));
verifyEqual(testCase,string(trial.suggest_categorical( ...
    "kind",["a","b"])),string(expected.values.kind));
verifyEqual(testCase,trial.Params.x,double(expected.values.x));
trial.report(2.5,3);
verifyEqual(testCase,trial.should_prune(),logical(expected.should_prune));
trial.set_user_attr("owner","matlab");
trial.set_system_attr("generation",2);
verifyEqual(testCase,string(trial.UserAttrs.owner), ...
    string(expected.user_attrs.owner));
verifyEqual(testCase,trial.SystemAttrs.generation, ...
    double(expected.system_attrs.generation));

lastwarn("");
repeated=trial.suggest_float("x",0.6,1);
[~,warningId]=lastwarn;
verifyEqual(testCase,repeated,double(expected.repeated_out_of_range));
verifyEqual(testCase,string(warningId),"radia:optuna:FixedParameter");
verifyEqual(testCase,string(trial.Distributions.x.name), ...
    string(expected.distribution_types.x));
verifyEqual(testCase,string(trial.Distributions.n.name), ...
    string(expected.distribution_types.n));
verifyEqual(testCase,string(trial.Distributions.kind.name), ...
    string(expected.distribution_types.kind));
verifyError(testCase,@()trial.suggest_float("missing",0,1), ...
    "radia:optuna:FixedParameterMissing");
verifyError(testCase,@()trial.suggest_int("x",0,2), ...
    "radia:optuna:IncompatibleDistribution");
invalid=radia.optuna.FixedTrial(struct("choice","z"));
verifyError(testCase,@()invalid.suggest_categorical( ...
    "choice",["a","b"]),"radia:optuna:FixedParameter");
end

function testFrozenTrialPublicBehaviorMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.frozen_trial;
trial=radia.optuna.create_trial(value=1.2, ...
    params=struct("x",0.5,"n",3,"kind","b"), ...
    distributions=struct( ...
        "x",radia.optuna.FloatDistribution(0,1), ...
        "n",radia.optuna.IntDistribution(1,5), ...
        "kind",radia.optuna.CategoricalDistribution(["a","b"])), ...
    intermediate_values=table(2,3,datetime("now",TimeZone="local"), ...
        VariableNames=["Step","Value","Timestamp"]), ...
    user_attrs=struct("owner","upstream"), ...
    system_attrs=struct("generation",1));
verifyEqual(testCase,trial.suggest_float("x",0,1),double(expected.values.x));
verifyEqual(testCase,trial.suggest_int("n",1,5),double(expected.values.n));
verifyEqual(testCase,string(trial.suggest_categorical( ...
    "kind",["a","b"])),string(expected.values.kind));
trial.report(99,9);
verifyEqual(testCase,trial.last_step(),double(expected.last_step));
verifyEqual(testCase,trial.should_prune(),logical(expected.should_prune));
trial.set_user_attr("owner","matlab");
trial.set_system_attr("generation",2);
verifyEqual(testCase,string(trial.UserAttrs.owner),string(expected.user_owner));
verifyEqual(testCase,trial.SystemAttrs.generation, ...
    double(expected.system_generation));
verifyError(testCase,@()trial.suggest_float("missing",0,1), ...
    "radia:optuna:FrozenParameterMissing");
verifyError(testCase,@()trial.suggest_int("x",0,2), ...
    "radia:optuna:IncompatibleDistribution");
verifyError(testCase,@()trial.suggest_categorical("kind",["a","c"]), ...
    "radia:optuna:IncompatibleDistribution");
end

function testStudyManagementMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.study_management;
source=string(tempname("C:\temp"))+".mat";
target=string(tempname("C:\temp"))+".mat";
multiTarget=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()cleanupStudyFiles([source,target,multiTarget]));
study=radia.optuna.create_study(study_name="alpha", ...
    direction="maximize",storage=source);
study.add_trial(radia.optuna.createTrial(Value=1.25, ...
    Params=struct("x",0.5), ...
    Distributions=struct("x",radia.optuna.FloatDistribution(0,1)), ...
    UserAttrs=struct("origin","oracle")));

loaded=radia.optuna.load_study(study_name="alpha",storage=source);
verifyEqual(testCase,loaded.Name,string(expected.loaded.name));
verifyEqual(testCase,string(loaded.direction()), ...
    string(expected.loaded.direction));
verifyEqual(testCase,height(loaded.TrialTable), ...
    double(expected.loaded.trial_count));
verifyEqual(testCase,loaded.best_value(),double(expected.loaded.best_value));
verifyEqual(testCase,radia.optuna.get_all_study_names(source), ...
    string(expected.names_before_delete));

radia.optuna.copy_study(from_study_name="alpha", ...
    from_storage=source,to_storage=target,to_study_name="beta");
copied=radia.optuna.load_study(study_name="beta",storage=target);
verifyEqual(testCase,copied.Name,string(expected.copied.name));
verifyEqual(testCase,string(copied.direction()), ...
    string(expected.copied.direction));
verifyEqual(testCase,copied.best_value(),double(expected.copied.best_value));
snapshot=copied.best_trial();
verifyEqual(testCase,snapshot.Params.x,double(expected.copied.param_x));
verifyEqual(testCase,string(snapshot.UserAttrs.origin), ...
    string(expected.copied.user_origin));

summary=radia.optuna.get_all_study_summaries(target);
verifyClass(testCase,summary,"radia.optuna.StudySummary");
verifyEqual(testCase,string(class(summary)), ...
    "radia.optuna."+string(expected.summary.type));
verifyEqual(testCase,string(summary.study_name), ...
    string(expected.summary.name));
verifyEqual(testCase,string(summary.directions), ...
    string(expected.summary.directions));
verifyEqual(testCase,string(summary.direction), ...
    string(expected.summary.direction));
verifyEqual(testCase,summary.n_trials,double(expected.summary.trial_count));
verifyEqual(testCase,summary.best_trial.Value, ...
    double(expected.summary.best_value));
withoutBest=radia.optuna.get_all_study_summaries( ...
    target,include_best_trial=false);
verifyEqual(testCase,~isempty(withoutBest.best_trial), ...
    logical(expected.summary.without_best_has_best));

radia.optuna.create_study(study_name="multi", ...
    directions=["minimize","maximize"],storage=multiTarget);
multiSummary=radia.optuna.get_all_study_summaries(multiTarget);
verifyEqual(testCase,string(multiSummary.directions), ...
    reshape(string(expected.multi_summary.directions),1,[]));
verifyEqual(testCase,string(expected.multi_summary.direction_error), ...
    "RuntimeError");
verifyError(testCase,@()readSummaryDirection(multiSummary), ...
    "radia:optuna:MultiObjectiveDirection");

radia.optuna.delete_study(study_name="alpha",storage=source);
verifyEqual(testCase,radia.optuna.get_all_study_names(source), ...
    strings(0,1));
verifyEqual(testCase,numel(expected.names_after_delete),0);
clear cleanup
cleanupStudyFiles([source,target,multiTarget]);
end

function testParameterImportancesMatchUpstream(testCase)
expected=testCase.TestData.Oracle.importance;
study=radia.optuna.Study(AutoSave=false);
for index=1:numel(expected.trials)
    row=expected.trials(index);
    study.add_trial(radia.optuna.createTrial(Value=double(row.value), ...
        Params=struct("x",double(row.x),"y",double(row.y), ...
        "mode",string(row.mode)), ...
        Distributions=struct( ...
        "x",radia.optuna.FloatDistribution(-1,1), ...
        "y",radia.optuna.FloatDistribution(0,1), ...
        "mode",radia.optuna.CategoricalDistribution(["A","B","C"]))));
end
result=radia.optuna.get_param_importances(study, ...
    evaluator=string(expected.evaluator.name), ...
    n_trees=double(expected.evaluator.n_trees), ...
    max_depth=double(expected.evaluator.max_depth), ...
    seed=double(expected.evaluator.seed));
verifyEqual(testCase,result.Parameter,string(expected.parameter_order));
verifyEqual(testCase,result.Importance, ...
    reshape(double(expected.values),[],1),AbsTol=0);
verifyEqual(testCase,sum(result.Importance),1,AbsTol=10*eps);
end

function testTerminationContractsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.terminator;
minimize=radia.optuna.Study(AutoSave=false);
for value=[5,4,4.5,4.2]
    trial=minimize.ask(); minimize.tell(trial,value);
end
maximize=radia.optuna.Study(Directions="maximize",AutoSave=false);
for value=[1,3,2.5,2]
    trial=maximize.ask(); maximize.tell(trial,value);
end
evaluator=radia.optuna.BestValueStagnationEvaluator(3);
verifyEqual(testCase,evaluator.evaluate( ...
    minimize.TrialTable,"minimize"),double(expected.remaining_minimize));
verifyEqual(testCase,evaluator.evaluate( ...
    maximize.TrialTable,"maximize"),double(expected.remaining_maximize));

maxTrials=radia.optuna.MaxTrialsCallback(3);
callbackStudy=radia.optuna.Study(AutoSave=false);
callbackStudy.optimize(@(trial)trial.Number,10, ...
    Callbacks=maxTrials.callback());
verifyEqual(testCase,height(callbackStudy.TrialTable), ...
    double(expected.max_trials_callback_count));

terminator=radia.optuna.Terminator( ...
    ImprovementEvaluator= ...
    radia.optuna.BestValueStagnationEvaluator(2),MinNTrials=1);
terminatorCallback=radia.optuna.TerminatorCallback(terminator);
values=reshape(double(expected.terminator_values),1,[]);
terminated=radia.optuna.Study(AutoSave=false);
terminated.optimize(@(trial)values(trial.Number+1),6, ...
    Callbacks=terminatorCallback.callback());
verifyEqual(testCase,height(terminated.TrialTable), ...
    double(expected.terminator_trial_count));
verifyEqual(testCase,terminated.TrialTable.Value,values');
end

function testRandomSamplerSeededSequence(testCase)
expected=testCase.TestData.Oracle.random_sampler_seed_123;
study=radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(123),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    actual=struct( ...
        "x",trial.suggest_float("x",-1,1), ...
        "q",trial.suggest_float("q",0,1,Step=0.2), ...
        "mesh",trial.suggest_int("mesh",1,9,Step=2), ...
        "log_mesh",trial.suggest_int("log_mesh",1,100,Log=true), ...
        "mode",trial.suggest_categorical("mode",["A","B","C"]));
    study.tell(trial,actual.x);
    verifyEqual(testCase,actual.x,expected(index).x,AbsTol=0);
    verifyEqual(testCase,actual.q,expected(index).q,AbsTol=0);
    verifyEqual(testCase,actual.mesh,double(expected(index).mesh));
    verifyEqual(testCase,actual.log_mesh,double(expected(index).log_mesh));
    verifyEqual(testCase,string(actual.mode),string(expected(index).mode));
end
end

function testDistributionCompatibilityMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.distributions;
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(4),AutoSave=false);
trial=study.ask();
warning("off","radia:optuna:DistributionAdjusted");
cleanup=onCleanup(@()warning("on","radia:optuna:DistributionAdjusted"));
value=trial.suggest_int("mesh order",2,10,Step=3);
lastwarn("");
repeated=trial.suggest_int("mesh order",2,10,Step=2);
[~,warningId]=lastwarn;
verifyEqual(testCase,value,double(expected.integer_value));
verifyEqual(testCase,repeated,double(expected.inconsistent_repeat_value));
verifyEqual(testCase,string(warningId),"radia:optuna:InconsistentParameter");
verifyError(testCase,@()trial.suggest_float("mesh order",2,10), ...
    "radia:optuna:IncompatibleDistribution");
verifyError(testCase,@()trial.suggest_int("mesh order",2,10,Log=true), ...
    "radia:optuna:IncompatibleDistribution");

categorical=study.ask();
categorical.suggest_categorical("choice",["A","B"]);
verifyError(testCase,@()categorical.suggest_categorical( ...
    "choice",["A","C"]),"radia:optuna:IncompatibleDistribution");

collision=study.ask();
first=collision.suggest_float("a-b",0,1);
second=collision.suggest_float("a_b",0,1);
verifyEqual(testCase,first,double(expected.colliding_names.a_b),AbsTol=0);
verifyEqual(testCase,second,double(expected.colliding_names.a_b_1),AbsTol=0);
verifyEqual(testCase,study.ParamTable.Name(end-1:end),["a-b";"a_b"]);

stepped=study.ask();
steppedValue=stepped.suggest_float("q",0,1,Step=0.3);
verifyEqual(testCase,steppedValue,double(expected.stepped_float.value),AbsTol=0);
verifyEqual(testCase,stepped.Distributions.q.high, ...
    double(expected.stepped_float.effective_high),AbsTol=2*eps);
clear cleanup
end

function testIntersectionSearchSpaceMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.search_space;
base=struct( ...
    "x",radia.optuna.FloatDistribution(0,1), ...
    "fixed",radia.optuna.FloatDistribution(2,2), ...
    "cat",radia.optuna.CategoricalDistribution(["A","B"]));
complete1=radia.optuna.create_trial(state="COMPLETE",value=1, ...
    params=struct("x",0.5,"fixed",2,"cat","A"), ...
    distributions=base);
complete2=radia.optuna.create_trial(state="COMPLETE",value=1, ...
    params=struct("x",0.5,"fixed",2,"cat","A","z",1), ...
    distributions=struct("x",base.x,"fixed",base.fixed, ...
    "cat",base.cat,"z",radia.optuna.IntDistribution(1,3)));
pruned=radia.optuna.create_trial(state="PRUNED", ...
    params=struct("x",0.5,"fixed",2,"cat","A"), ...
    distributions=struct( ...
    "x",radia.optuna.FloatDistribution(-1,1), ...
    "fixed",base.fixed,"cat",base.cat));
failed=radia.optuna.create_trial(state="FAIL", ...
    params=struct("x",0.5),distributions=struct( ...
    "x",radia.optuna.FloatDistribution(-2,2)));
waiting=radia.optuna.FrozenTrial(State="WAITING", ...
    Params=struct("x",0.5),Distributions=struct( ...
    "x",radia.optuna.FloatDistribution(-3,3)));
trials=[complete1;complete2;pruned;failed;waiting];

withoutPruned=radia.optuna.intersection_search_space(trials);
withPruned=radia.optuna.intersection_search_space( ...
    trials,include_pruned=true);
verifyEqual(testCase,sort(string(keys(withoutPruned))), ...
    sort(reshape(string(expected.without_pruned),1,[])));
verifyEqual(testCase,sort(string(keys(withPruned))), ...
    sort(reshape(string(expected.with_pruned),1,[])));
verifyTrue(testCase,logical(expected.single_distribution_is_included));
verifyTrue(testCase,isKey(withoutPruned,"fixed"));

study=radia.optuna.Study(AutoSave=false);
study.add_trials(trials);
calculator=radia.optuna.IntersectionSearchSpace();
calculated=calculator.calculate(study);
verifyEqual(testCase,sort(string(keys(calculated))), ...
    sort(reshape(string(expected.calculator),1,[])));
end

function testStudyDirectionAndTrialStateMatchUpstream(testCase)
expected=testCase.TestData.Oracle.enums;
directions=[radia.optuna.StudyDirection.NOT_SET, ...
    radia.optuna.StudyDirection.MINIMIZE, ...
    radia.optuna.StudyDirection.MAXIMIZE];
states=[radia.optuna.TrialState.RUNNING, ...
    radia.optuna.TrialState.COMPLETE,radia.optuna.TrialState.PRUNED, ...
    radia.optuna.TrialState.FAIL,radia.optuna.TrialState.WAITING];
verifyEqual(testCase,string(directions), ...
    reshape(string({expected.study_direction.name}),1,[]));
verifyEqual(testCase,double(directions), ...
    reshape(double([expected.study_direction.value]),1,[]));
verifyEqual(testCase,string(states), ...
    reshape(string({expected.trial_state.name}),1,[]));
verifyEqual(testCase,double(states), ...
    reshape(double([expected.trial_state.value]),1,[]));
verifyEqual(testCase,states.is_finished(), ...
    reshape(logical([expected.trial_state.is_finished]),1,[]));

study=radia.optuna.create_study( ...
    direction=radia.optuna.StudyDirection.MINIMIZE,AutoSave=false);
verifyEqual(testCase,study.direction(), ...
    radia.optuna.StudyDirection.MINIMIZE);
trial=radia.optuna.create_trial( ...
    state=radia.optuna.TrialState.PRUNED,value=3);
verifyEqual(testCase,trial.State,"PRUNED");
end

function testUnfinishedAddedTrialsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.unfinished_trials;
distribution=struct("x",radia.optuna.FloatDistribution(0,1));
waiting=radia.optuna.create_trial(state=radia.optuna.TrialState.WAITING, ...
    params=struct("x",0.25),distributions=distribution, ...
    user_attrs=struct("source","waiting"));
running=radia.optuna.create_trial(state=radia.optuna.TrialState.RUNNING, ...
    params=struct("x",0.75),distributions=distribution, ...
    user_attrs=struct("source","running"));

verifyEqual(testCase,waiting.Number,double(expected.factory.waiting_number));
verifyEqual(testCase,~isnat(waiting.DatetimeStart), ...
    logical(expected.factory.waiting_has_start));
verifyEqual(testCase,~isnat(waiting.DatetimeComplete), ...
    logical(expected.factory.waiting_has_complete));
verifyEqual(testCase,~isnan(waiting.Duration), ...
    logical(expected.factory.waiting_has_duration));
verifyEqual(testCase,running.Number,double(expected.factory.running_number));
verifyEqual(testCase,~isnat(running.DatetimeStart), ...
    logical(expected.factory.running_has_start));
verifyEqual(testCase,~isnat(running.DatetimeComplete), ...
    logical(expected.factory.running_has_complete));
verifyEqual(testCase,~isnan(running.Duration), ...
    logical(expected.factory.running_has_duration));

study=radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(9),AutoSave=false);
study.add_trials([waiting;running]);
before=study.get_trials();
verifyEqual(testCase,[before.Number], ...
    reshape(double([expected.before.number]),1,[]));
verifyEqual(testCase,[before.State], ...
    reshape(string({expected.before.state}),1,[]));
verifyEqual(testCase,[before.Params],struct("x",{0.25,0.75}));
verifyEqual(testCase,~isnat([before.DatetimeStart]), ...
    reshape(logical([expected.before.has_start]),1,[]));
verifyEqual(testCase,~isnat([before.DatetimeComplete]), ...
    reshape(logical([expected.before.has_complete]),1,[]));
verifyEqual(testCase,~isnan([before.Duration]), ...
    reshape(logical([expected.before.has_duration]),1,[]));

claimed=study.ask();
verifyEqual(testCase,claimed.Number,double(expected.claimed_number));
verifyEqual(testCase,claimed.suggest_float("x",0,1), ...
    double(expected.claimed_value),AbsTol=0);
verifyEqual(testCase,claimed.Params.x,double(expected.claimed_params.x));
verifyEqual(testCase,string(claimed.UserAttrs.source), ...
    string(expected.claimed_user_attrs.source));
afterClaim=study.get_trials();
verifyEqual(testCase,[afterClaim.State], ...
    reshape(string(expected.after_claim_states),1,[]));
verifyEqual(testCase,~isnat([afterClaim.DatetimeStart]), ...
    reshape(logical(expected.after_claim_has_start),1,[]));

fresh=study.ask();
verifyEqual(testCase,fresh.Number,double(expected.fresh_number));
verifyTrue(testCase,isempty(fieldnames(fresh.Params)));
end

function testBaseClassHierarchyMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.base_trial;
fixed=radia.optuna.FixedTrial(struct("x",0.5));
frozen=radia.optuna.create_trial(value=1);
study=radia.optuna.Study(AutoSave=false);
trial=study.ask();
actual=[isa(fixed,"radia.optuna.BaseTrial"), ...
    isa(frozen,"radia.optuna.BaseTrial"), ...
    isa(trial,"radia.optuna.BaseTrial")];
verifyEqual(testCase,actual,[logical(expected.is_base_trial.fixed), ...
    logical(expected.is_base_trial.frozen), ...
    logical(expected.is_base_trial.trial)]);
verifyEqual(testCase,[fixed.number(),frozen.number(),trial.number()], ...
    [double(expected.numbers.fixed),double(expected.numbers.frozen), ...
    double(expected.numbers.trial)]);
verifyEqual(testCase,string(expected.construction_error),"TypeError");
verifyError(testCase,@()radia.optuna.BaseTrial(), ...
    "MATLAB:class:abstract");

components=testCase.TestData.Oracle.base_components;
sampler=radia.optuna.RandomSampler(5);
pruner=radia.optuna.NopPruner();
verifyEqual(testCase,isa(sampler,"radia.optuna.BaseSampler"), ...
    logical(components.sampler_is_base));
verifyEqual(testCase,isa(pruner,"radia.optuna.BasePruner"), ...
    logical(components.pruner_is_base));
verifyEqual(testCase,pruner.prune(study,trial), ...
    logical(components.nop_decision));
verifyEqual(testCase,string(components.construction_errors.sampler), ...
    "TypeError");
verifyEqual(testCase,string(components.construction_errors.pruner), ...
    "TypeError");
verifyError(testCase,@()radia.optuna.BaseSampler(), ...
    "MATLAB:class:abstract");
verifyError(testCase,@()radia.optuna.BasePruner(), ...
    "MATLAB:class:abstract");
end

function testDistributionJSONAndCompatibilityMatchUpstream(testCase)
expected=testCase.TestData.Oracle.distribution_json;
warning("off","radia:optuna:FutureWarning");
cleanup=onCleanup(@()warning("on","radia:optuna:FutureWarning"));
actual=struct( ...
    "float",radia.optuna.FloatDistribution(0,1), ...
    "log_float",radia.optuna.FloatDistribution(0.001,1,Log=true), ...
    "stepped_float",radia.optuna.FloatDistribution(0,1,Step=0.2), ...
    "integer",radia.optuna.IntDistribution(1,9,Step=2), ...
    "log_integer",radia.optuna.IntDistribution(1,100,Log=true), ...
    "categorical",radia.optuna.CategoricalDistribution({"A",2,true}), ...
    "uniform",radia.optuna.UniformDistribution(0,1), ...
    "log_uniform",radia.optuna.LogUniformDistribution(0.001,1), ...
    "discrete_uniform",radia.optuna.DiscreteUniformDistribution(0,1,0.2), ...
    "int_uniform",radia.optuna.IntUniformDistribution(1,9,2), ...
    "int_log_uniform",radia.optuna.IntLogUniformDistribution(1,100,1));
names=string(fieldnames(actual));
for index=1:numel(names)
    name=names(index);
    encoded=radia.optuna.distribution_to_json(actual.(name));
    verifyEqual(testCase,jsondecode(encoded), ...
        jsondecode(expected.encoded.(name)));
    decoded=radia.optuna.json_to_distribution(encoded);
    verifyEqual(testCase,string(decoded.name), ...
        string(expected.roundtrip_types.(name)));
end
verifyTrue(testCase,radia.optuna.check_distribution_compatibility( ...
    radia.optuna.FloatDistribution(0,1), ...
    radia.optuna.FloatDistribution(-1,2,Step=0.2)));
verifyTrue(testCase,logical(expected.compatibility.range_change_allowed));
verifyEqual(testCase,string(struct2cell(expected.compatibility.errors)), ...
    repmat("ValueError",3,1));
verifyError(testCase,@()radia.optuna.check_distribution_compatibility( ...
    radia.optuna.FloatDistribution(0,1), ...
    radia.optuna.IntDistribution(0,1)), ...
    "radia:optuna:IncompatibleDistribution");
verifyError(testCase,@()radia.optuna.check_distribution_compatibility( ...
    radia.optuna.FloatDistribution(0,1), ...
    radia.optuna.FloatDistribution(0.001,1,Log=true)), ...
    "radia:optuna:IncompatibleDistribution");
verifyError(testCase,@()radia.optuna.check_distribution_compatibility( ...
    radia.optuna.CategoricalDistribution(["A","B"]), ...
    radia.optuna.CategoricalDistribution(["A","C"])), ...
    "radia:optuna:IncompatibleDistribution");
clear cleanup
end

function testTPESamplerSeededSequence(testCase)
expected=reshape(double(testCase.TestData.Oracle.tpe_sampler_seed_37),[],1);
study=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=37,NStartupTrials=4),AutoSave=false);
actual=zeros(size(expected));
for index=1:numel(expected)
    trial=study.ask();
    actual(index)=trial.suggest_float("x",-2,2);
    study.tell(trial,(actual(index)-0.25)^2);
end
verifyEqual(testCase,actual,expected,AbsTol=5e-12);
end

function testTPECallableGammaWeightsSeededSequence(testCase)
expected=testCase.TestData.Oracle.custom_tpe_sampler_gamma_weights;
sampler=radia.optuna.TPESampler(Seed=97,NStartupTrials=4, ...
    GammaFcn=@(count)min(3,count),WeightsFcn=@customTPEWeights);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
single=reshape(double(expected.single),[],1);
actual=zeros(size(single));
for index=1:numel(single)
    trial=study.ask();
    actual(index)=trial.suggest_float("x",-2,2);
    study.tell(trial,(actual(index)-0.35)^2);
end
verifyEqual(testCase,actual,single,AbsTol=5e-12);

sampler=radia.optuna.TPESampler(Seed=101,NStartupTrials=4, ...
    GammaFcn=@(count)min(3,count),WeightsFcn=@customTPEWeights);
study=radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
multi=expected.multi;
for index=1:numel(multi)
    trial=study.ask();
    x=trial.suggest_float("x",-2,2);
    y=trial.suggest_float("y",-1,3);
    study.tell(trial,[(x-0.4)^2+0.1*y*y,(y+0.2)^2]);
    verifyEqual(testCase,x,double(multi(index).x),AbsTol=5e-12);
    verifyEqual(testCase,y,double(multi(index).y),AbsTol=5e-12);
end
end

function testTPEGroupAndIndependentWarningMatchUpstream(testCase)
contract=testCase.TestData.Oracle.tpe_group;
expected=contract.sequence;
study=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=101,NStartupTrials=4,Multivariate=true,Group=true), ...
    AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    branch=string(trial.suggest_categorical( ...
        "branch",["left","right"]));
    x=trial.suggest_float("x",-1,1);
    if branch=="left"
        y=trial.suggest_float("y",0,2);
        value=(x-0.2)^2+(y-0.4)^2;
        verifyEqual(testCase,y,double(expected(index).y),AbsTol=5e-12);
    else
        z=trial.suggest_int("z",1,5);
        value=(x+0.1)^2+0.05*z;
        verifyEqual(testCase,z,double(expected(index).z));
    end
    study.tell(trial,value);
    verifyEqual(testCase,trial.Number,double(expected(index).number));
    verifyEqual(testCase,branch,string(expected(index).branch));
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
end

verifyGreaterThan(testCase, ...
    double(contract.independent_warning_enabled_count),0);
verifyEqual(testCase, ...
    double(contract.independent_warning_disabled_count),0);
verifyEqual(testCase,double(contract.group_warning_count),0);

enabled=tpeWarningStudy(true,false);
trial=enabled.ask();
verifyWarning(testCase,@()trial.suggest_float("y",0,1), ...
    "radia:optuna:TPEIndependentSampling");

disabled=tpeWarningStudy(false,false);
trial=disabled.ask();
verifyWarningFree(testCase,@()trial.suggest_float("y",0,1));

grouped=tpeWarningStudy(true,true);
trial=grouped.ask();
verifyWarningFree(testCase,@()trial.suggest_float("y",0,1));
end

function testTPECategoricalDistanceMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.tpe_categorical_distance;
distances=containers.Map('KeyType','char','ValueType','any');
distances('level')=@categoricalLevelDistance;
sampler=radia.optuna.TPESampler(Seed=107,NStartupTrials=4, ...
    CategoricalDistanceFcn=distances);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
levels=["zero","one","two","three"];
for index=1:numel(expected)
    trial=study.ask();
    level=string(trial.suggest_categorical("level",levels));
    position=find(levels==level,1)-1;
    study.tell(trial,(position-1.3)^2);
    verifyEqual(testCase,trial.Number,double(expected(index).number));
    verifyEqual(testCase,level,string(expected(index).level));
    verifyEqual(testCase,position,double(expected(index).position));
end
end

function testMultiObjectiveTPESamplerSeededSequence(testCase)
expected=testCase.TestData.Oracle.multiobjective_tpe_sampler_seed_41;
study=radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.TPESampler(Seed=41,NStartupTrials=4), ...
    AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-2,2);
    y=trial.suggest_float("y",-1,3);
    study.tell(trial,[(x-0.4)^2+0.1*y*y, ...
        (y+0.2)^2+0.1*x*x]);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
    verifyEqual(testCase,y,double(expected(index).y),AbsTol=5e-12);
end
end

function study=tpeWarningStudy(warn,group)
sampler=radia.optuna.TPESampler(Seed=103,NStartupTrials=1, ...
    Multivariate=true,Group=group,WarnIndependentSampling=warn);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
first=study.ask();
x=first.suggest_float("x",0,1);
study.tell(first,x);
second=study.ask();
x=second.suggest_float("x",0,1);
y=second.suggest_float("y",0,1);
study.tell(second,x+y);
end

function distance=categoricalLevelDistance(first,second)
levels=["zero","one","two","three"];
distance=abs(find(levels==string(first),1)- ...
    find(levels==string(second),1));
end

function testMixedTPESamplerSeededSequence(testCase)
expected=testCase.TestData.Oracle.mixed_tpe_sampler_seed_43;
study=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=43,NStartupTrials=4),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    q=trial.suggest_float("q",0,1,Step=0.2);
    mesh=trial.suggest_int("mesh",1,9,Step=2);
    logMesh=trial.suggest_int("log_mesh",1,100,Log=true);
    mode=trial.suggest_categorical("mode",["A","B","C"]);
    modePenalty=0.2*double(mode=="A")+0.35*double(mode=="C");
    loss=(x-0.2)^2+(q-0.6)^2+0.01*mesh+0.001*logMesh+modePenalty;
    study.tell(trial,loss);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
    verifyEqual(testCase,q,double(expected(index).q),AbsTol=5e-12);
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,logMesh,double(expected(index).log_mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testGridSamplerSeededOrder(testCase)
expected=testCase.TestData.Oracle.grid_sampler_seed_17;
space=struct("x",[-1,0,1],"mode",["A","B"]);
study=radia.optuna.Study( ...
    Sampler=radia.optuna.GridSampler(space,Seed=17),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,x);
    verifyEqual(testCase,x,double(expected(index).x));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testNSGAIISeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.nsgaii_sampler_seed_19;
study=radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIISampler(Seed=19,PopulationSize=4), ...
    AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mesh=trial.suggest_int("mesh",1,5,Step=2);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,[x*x+0.1*mesh, ...
        (x-0.5)^2+0.2*double(mode=="B")]);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testNSGAIIISeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.nsgaiii_sampler_seed_23;
study=radia.optuna.Study( ...
    Directions=["minimize","minimize","minimize"], ...
    Sampler=radia.optuna.NSGAIIISampler(Seed=23,PopulationSize=4), ...
    AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mesh=trial.suggest_int("mesh",1,5,Step=2);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,[x*x+0.1*mesh, ...
        (x-0.5)^2+0.2*double(mode=="B"), ...
        (x+0.25)^2+0.05*mesh]);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testBruteForceSeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.brute_force_sampler_seed_29;
study=radia.optuna.Study( ...
    Sampler=radia.optuna.BruteForceSampler(Seed=29),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    mesh=trial.suggest_int("mesh",1,3);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,mesh+0.1*double(mode=="B"));
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testConditionalBruteForceSeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.conditional_brute_force_sampler_seed_79;
study=radia.optuna.Study(Sampler= ...
    radia.optuna.BruteForceSampler(Seed=79),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    branch=trial.suggest_categorical("branch",["depth","mode"]);
    if branch=="depth"
        value=string(trial.suggest_int("depth",1,2));
        parameter="depth";
        loss=double(value);
    else
        value=string(trial.suggest_categorical("mode",["A","B","C"]));
        parameter="mode";
        loss=double(char(value)-'A');
    end
    study.tell(trial,loss);
    verifyEqual(testCase,branch,string(expected(index).branch));
    verifyEqual(testCase,parameter,string(expected(index).parameter));
    verifyEqual(testCase,value,string(expected(index).value));
end
end

function testCmaEsSeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.cmaes_sampler_seed_31;
study=radia.optuna.Study(Sampler=radia.optuna.CmaEsSampler( ...
    Seed=31,NStartupTrials=1,PopulationSize=4),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-2,2);
    y=trial.suggest_float("y",-1,3);
    study.tell(trial,(x-0.4)^2+0.5*(y+0.2)^2);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
    verifyEqual(testCase,y,double(expected(index).y),AbsTol=5e-12);
end
end

function testCmaEsIndependentSamplerSeededSequence(testCase)
expected=testCase.TestData.Oracle.cmaes_independent_sampler_seed_31;
sampler=radia.optuna.CmaEsSampler(Seed=31,NStartupTrials=1, ...
    PopulationSize=4,IndependentSampler=radia.optuna.RandomSampler(211), ...
    WarnIndependentSampling=false);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-2,2);
    y=trial.suggest_float("y",-1,3);
    mode=string(trial.suggest_categorical("mode",["A","B","C"]));
    penalty=0.1*find(["A","B","C"]==mode,1)-0.1;
    study.tell(trial,(x-0.4)^2+0.5*(y+0.2)^2+penalty);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
    verifyEqual(testCase,y,double(expected(index).y),AbsTol=5e-12);
    verifyEqual(testCase,mode,string(expected(index).mode));
end
end

function testScrambledQMCSeededProposalSequence(testCase)
contract=testCase.TestData.Oracle.scrambled_qmc_sampler_seed_47;
for type=["sobol","halton"]
    expected=contract.(type);
    study=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
        QMCType=type,Scramble=true,Seed=47),AutoSave=false);
    for index=1:numel(expected)
        trial=study.ask();
        x=trial.suggest_float("x",-1,1);
        y=trial.suggest_float("y",0,4);
        study.tell(trial,x*x+y*y);
        verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
        verifyEqual(testCase,y,double(expected(index).y),AbsTol=0);
    end
end
end

function testGPSeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.gp_sampler_seed_53;
study=radia.optuna.Study(Sampler=radia.optuna.GPSampler( ...
    Seed=53,NStartupTrials=10),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mesh=trial.suggest_int("mesh",1,5,Step=2);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,x*x+0.1*mesh+0.2*double(mode=="B"));
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testGPConstraintsSeededProposalSequence(testCase)
expected=testCase.TestData.Oracle.gp_constraints_sampler_seed_89;
sampler=radia.optuna.GPSampler(Seed=89,NStartupTrials=5, ...
    ConstraintsFcn=@(trial)trial.UserAttrs.constraints);
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    trial.set_user_attr("constraints",x-0.1);
    study.tell(trial,(x-0.35)^2);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
    verifyEqual(testCase,study.constraintsForTrial(trial.Number), ...
        double(expected(index).constraint),AbsTol=0);
end
end

function testGPTableResumeReplaysUpstreamHistory(testCase)
expected=testCase.TestData.Oracle.gp_sampler_seed_53;
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()cleanupStudyFiles(path));
study=radia.optuna.Study(StoragePath=path,Sampler= ...
    radia.optuna.GPSampler(Seed=53,NStartupTrials=10));
for index=1:numel(expected)-1
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mesh=trial.suggest_int("mesh",1,5,Step=2);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,x*x+0.1*mesh+0.2*double(mode=="B"));
end
clear study
resumed=radia.optuna.load_study(storage=path,sampler= ...
    radia.optuna.GPSampler(Seed=53,NStartupTrials=10));
trial=resumed.ask();
x=trial.suggest_float("x",-1,1);
mesh=trial.suggest_int("mesh",1,5,Step=2);
mode=trial.suggest_categorical("mode",["A","B"]);
last=expected(end);
verifyEqual(testCase,x,double(last.x),AbsTol=0);
verifyEqual(testCase,mesh,double(last.mesh));
verifyEqual(testCase,string(mode),string(last.mode));
clear cleanup
cleanupStudyFiles(path);
end

function testNSGAIICrossoverSeededSequences(testCase)
expected=testCase.TestData.Oracle.nsgaii_crossovers_seed_73;
cases={ ...
    "uniform",radia.optuna.nsgaii.UniformCrossover(); ...
    "blxalpha",radia.optuna.nsgaii.BLXAlphaCrossover(); ...
    "sbx",radia.optuna.nsgaii.SBXCrossover(); ...
    "vsbx",radia.optuna.nsgaii.VSBXCrossover(); ...
    "spx",radia.optuna.nsgaii.SPXCrossover(); ...
    "undx",radia.optuna.nsgaii.UNDXCrossover()};
for caseIndex=1:size(cases,1)
    name=cases{caseIndex,1};
    rows=expected.(name);
    sampler=radia.optuna.NSGAIISampler(Seed=73,PopulationSize=4, ...
        MutationProbability=0,CrossoverProbability=1, ...
        Crossover=cases{caseIndex,2});
    study=radia.optuna.Study(Sampler=sampler, ...
        Directions=["minimize","minimize"],AutoSave=false);
    for index=1:numel(rows)
        trial=study.ask();
        x=trial.suggest_float("x",-1,1);
        y=trial.suggest_float("y",-2,2);
        z=trial.suggest_float("z",0,3);
        study.tell(trial,[x*x+0.2*y*y+0.1*z, ...
            (x-0.4)^2+(y+0.3)^2+z*z]);
        verifyEqual(testCase,x,double(rows(index).x), ...
            sprintf("%s trial %d x",name,index),AbsTol=5e-12);
        verifyEqual(testCase,y,double(rows(index).y), ...
            sprintf("%s trial %d y",name,index),AbsTol=5e-12);
        verifyEqual(testCase,z,double(rows(index).z), ...
            sprintf("%s trial %d z",name,index),AbsTol=5e-12);
    end
end
end

function testPartialFixedSamplerSeededSequence(testCase)
expected=testCase.TestData.Oracle.partial_fixed_sampler_seed_61;
sampler=radia.optuna.PartialFixedSampler(struct("x",0.25), ...
    radia.optuna.RandomSampler(61));
study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",0,1);
    y=trial.suggest_float("y",-1,1);
    study.tell(trial,x*x+y*y);
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
    verifyEqual(testCase,y,double(expected(index).y),AbsTol=0);
end
end

function testMultivariateTPESeededSequence(testCase)
expected=testCase.TestData.Oracle.multivariate_tpe_sampler_seed_67;
study=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=67,NStartupTrials=4,Multivariate=true),AutoSave=false);
for index=1:numel(expected)
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    mesh=trial.suggest_int("mesh",1,5,Step=2);
    mode=trial.suggest_categorical("mode",["A","B"]);
    study.tell(trial,(x-0.2)^2+0.05*mesh+0.2*double(mode~="B"));
    verifyEqual(testCase,x,double(expected(index).x),AbsTol=5e-12);
    verifyEqual(testCase,mesh,double(expected(index).mesh));
    verifyEqual(testCase,string(mode),string(expected(index).mode));
end
end

function testUnscrambledQMCSeededSequence(testCase)
contract=testCase.TestData.Oracle.unscrambled_qmc_sampler_seed_71;
for type=["sobol","halton"]
    expected=contract.(type);
    study=radia.optuna.Study(Sampler=radia.optuna.QMCSampler( ...
        QMCType=type,Scramble=false,Seed=71),AutoSave=false);
    for index=1:numel(expected)
        trial=study.ask();
        x=trial.suggest_float("x",-1,1);
        y=trial.suggest_float("y",0,4);
        study.tell(trial,x*x+y*y);
        verifyEqual(testCase,x,double(expected(index).x),AbsTol=0);
        verifyEqual(testCase,y,double(expected(index).y),AbsTol=0);
    end
end
end

function testQMCWarningOptionsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.qmc_warnings;
verifyGreaterThan(testCase,double(expected.asynchronous_enabled_count),0);
verifyEqual(testCase,double(expected.asynchronous_disabled_count),0);
verifyGreaterThan(testCase,double(expected.independent_enabled_count),0);
verifyEqual(testCase,double(expected.independent_disabled_count),0);

verifyWarning(testCase,@()radia.optuna.QMCSampler( ...
    Scramble=true,WarnAsynchronousSeeding=true), ...
    "radia:optuna:QMCAsynchronousSeeding");
verifyWarningFree(testCase,@()radia.optuna.QMCSampler( ...
    Scramble=true,WarnAsynchronousSeeding=false));

enabled=radia.optuna.QMCSampler(Seed=11,WarnIndependentSampling=true);
enabledStudy=radia.optuna.Study(Sampler=enabled,AutoSave=false);
first=enabledStudy.ask();
first.suggest_categorical("kind",["a","b"]);
enabledStudy.tell(first,0);
second=enabledStudy.ask();
verifyWarning(testCase,@()second.suggest_categorical("kind",["a","b"]), ...
    "radia:optuna:QMCIndependentSampling");

disabled=radia.optuna.QMCSampler(Seed=11,WarnIndependentSampling=false);
disabledStudy=radia.optuna.Study(Sampler=disabled,AutoSave=false);
first=disabledStudy.ask();
first.suggest_categorical("kind",["a","b"]);
disabledStudy.tell(first,0);
second=disabledStudy.ask();
verifyWarningFree(testCase, ...
    @()second.suggest_categorical("kind",["a","b"]));
end

function testPrunerDecisionsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.pruners;
percentile=radia.optuna.Study(Pruner=radia.optuna.PercentilePruner( ...
    50,NStartupTrials=0,NWarmupSteps=0,IntervalSteps=2,NMinTrials=1), ...
    AutoSave=false);
addCompletedTrial(percentile,[1,3],[1,1]);
addCompletedTrial(percentile,[1,3],[3,3]);
trial=percentile.ask();
trial.report(5,1); trial.report(4,3);
verifyEqual(testCase,trial.shouldPrune(), ...
    logical(expected.percentile_minimize));

maximize=radia.optuna.Study(Directions="maximize",Pruner= ...
    radia.optuna.PercentilePruner(50,NStartupTrials=0),AutoSave=false);
addCompletedTrial(maximize,0,1); addCompletedTrial(maximize,0,3);
trial=maximize.ask(); trial.report(1.5,0);
verifyEqual(testCase,trial.shouldPrune(), ...
    logical(expected.percentile_maximize));

threshold=radia.optuna.Study(Pruner=radia.optuna.ThresholdPruner( ...
    Lower=0,Upper=10),AutoSave=false);
trial=threshold.ask(); trial.report(5,0);
actual(1)=trial.shouldPrune();
trial.report(11,1); actual(2)=trial.shouldPrune();
nanTrial=threshold.ask(); nanTrial.report(NaN,0);
actual(3)=nanTrial.shouldPrune();
verifyEqual(testCase,actual,reshape(logical(expected.threshold),1,[]));

patient=radia.optuna.Study(Pruner=radia.optuna.PatientPruner([], ...
    Patience=1,MinDelta=0),AutoSave=false);
trial=patient.ask();
patientValues=[10,4,5,6];
for index=1:4, trial.report(patientValues(index),index-1); end
wrapped=radia.optuna.Study(Pruner=radia.optuna.PatientPruner( ...
    radia.optuna.NopPruner(),Patience=1),AutoSave=false);
wrappedTrial=wrapped.ask();
for index=1:4, wrappedTrial.report(patientValues(index),index-1); end
verifyEqual(testCase,[trial.shouldPrune(),wrappedTrial.shouldPrune()], ...
    reshape(logical(expected.patient),1,[]));

halvingPruner=radia.optuna.SuccessiveHalvingPruner( ...
    MinResource=1,ReductionFactor=2);
halving=radia.optuna.Study(Pruner=halvingPruner,AutoSave=false);
first=halving.ask(); first.report(1,1);
firstDecision=first.shouldPrune(); halving.tell(first,1);
second=halving.ask(); second.report(2,1);
secondDecision=second.shouldPrune();
verifyEqual(testCase,[firstDecision,secondDecision], ...
    reshape(logical(expected.successive_halving.decisions),1,[]));
verifyEqual(testCase,[first.SystemAttrs.completed_rung_0, ...
    second.SystemAttrs.completed_rung_0], ...
    reshape(double(expected.successive_halving.rung_values),1,[]));
bootstrap=radia.optuna.Study(Pruner= ...
    radia.optuna.SuccessiveHalvingPruner(MinResource=1, ...
    ReductionFactor=2,BootstrapCount=1),AutoSave=false);
bootstrapTrial=bootstrap.ask(); bootstrapTrial.report(1,1);
verifyEqual(testCase,bootstrapTrial.shouldPrune(), ...
    logical(expected.successive_halving.bootstrap));

hyperbandPruner=radia.optuna.HyperbandPruner(MinResource=1, ...
    MaxResource=9,ReductionFactor=3);
hyperband=radia.optuna.Study(Name="hb",Pruner=hyperbandPruner,AutoSave=false);
hyperbandTrial=hyperband.ask(); hyperbandTrial.report(1,0);
verifyEqual(testCase,hyperbandTrial.shouldPrune(), ...
    logical(expected.hyperband.first_decision));
brackets=zeros(1,10);
for number=0:9, brackets(number+1)=hyperbandPruner.bracketId(hyperband,number); end
verifyEqual(testCase,brackets, ...
    reshape(double(expected.hyperband.bracket_ids),1,[]));

wilcoxon=radia.optuna.Study(Pruner=radia.optuna.WilcoxonPruner( ...
    PThreshold=0.1,NStartupSteps=2),AutoSave=false);
addCompletedTrial(wilcoxon,0:5,zeros(1,6));
wilcoxonTrial=wilcoxon.ask();
for step=0:5, wilcoxonTrial.report(10,step); end
nonfinite=wilcoxon.ask(); nonfinite.report(Inf,0);
warning("off","radia:optuna:WilcoxonNonfinite");
cleanup=onCleanup(@()warning("on","radia:optuna:WilcoxonNonfinite"));
verifyEqual(testCase,[wilcoxonTrial.shouldPrune(),nonfinite.shouldPrune()], ...
    reshape(logical(expected.wilcoxon),1,[]));
clear cleanup
end

function testConstraintParetoMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.constraints;
sampler=radia.optuna.NSGAIISampler(Seed=73,PopulationSize=4, ...
    ConstraintsFcn=@(trial)trial.UserAttrs.c);
study=radia.optuna.Study(Directions=["minimize","minimize"], ...
    Sampler=sampler,AutoSave=false);
values=[0,0;1,2;2,1;-1,-1];
constraints={1,-1,0,2};
for index=1:4
    trial=study.ask();
    trial.set_user_attr("c",constraints{index});
    study.tell(trial,values(index,:));
end
front=study.best_trials();
verifyEqual(testCase,sort(reshape([front.Number],[],1)), ...
    reshape(double(expected.pareto_trial_numbers),[],1));
verifyEqual(testCase,study.TrialTable.State, ...
    reshape(string(expected.states),[],1));
for index=1:4
    verifyEqual(testCase,study.constraintsForTrial(index-1), ...
        reshape(double(expected.constraints(index,:)),1,[]));
end
end

function addCompletedTrial(study,steps,values)
trial=study.ask();
for index=1:numel(steps), trial.report(values(index),steps(index)); end
study.tell(trial,values(end));
end

function value=failingConstraintCallback()
value=[]; %#ok<NASGU>
error("radia:test:ConstraintCallback","constraint callback failed");
end

function value=readSummaryDirection(summary)
value=summary.direction;
end

function cleanupStudyFiles(paths)
for path=reshape(string(paths),1,[])
    if isfile(path), delete(path); end
    if isfile(path+".bak"), delete(path+".bak"); end
end
end

function weights=customTPEWeights(count)
if count==0
    weights=zeros(0,1);
elseif count==1
    % numpy.linspace(0.2,1.0,1) returns the start value.
    weights=0.2;
else
    weights=linspace(0.2,1,count)';
end
end
