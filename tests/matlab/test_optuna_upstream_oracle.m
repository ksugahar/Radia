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
verifyEqual(testCase,radia.optuna.version(),string(oracle.optuna_version));
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

function testLoggingPublicContractMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.logging;
constants=["CRITICAL","DEBUG","ERROR","FATAL","INFO","WARN","WARNING"];
for name=constants
    actual=feval("radia.optuna."+name);
    verifyEqual(testCase,actual,double(expected.constants.(name)));
end
radia.optuna.set_verbosity(radia.optuna.INFO());
verifyEqual(testCase,radia.optuna.get_verbosity(),double(expected.initial));
logger=radia.optuna.get_logger("unit");
verifyEqual(testCase,logger.name,string(expected.logger.name));
verifyEqual(testCase,logger.level,double(expected.logger.level));
verifyEqual(testCase,logger.propagate,logical(expected.logger.propagate));
verifyEqual(testCase,logger.handlers,double(expected.logger.handlers));
formatter=radia.optuna.create_default_formatter();
verifyEqual(testCase,formatter.format,string(expected.formatter.format));
radia.optuna.set_verbosity(radia.optuna.DEBUG());
verifyEqual(testCase,radia.optuna.get_verbosity(),double(expected.after_set));
radia.optuna.disable_default_handler();
radia.optuna.disable_propagation();
root=radia.optuna.get_logger("optuna");
verifyEqual(testCase,root.propagate,logical(expected.disabled.propagate));
verifyEqual(testCase,root.handlers,double(expected.disabled.handlers));
radia.optuna.enable_default_handler();
radia.optuna.enable_propagation();
root=radia.optuna.get_logger("optuna");
verifyEqual(testCase,root.propagate,logical(expected.enabled.propagate));
verifyEqual(testCase,root.handlers,double(expected.enabled.handlers));
radia.optuna.set_verbosity(radia.optuna.WARNING());
radia.optuna.disable_propagation();
end

function testSamplerPublicMembersMatchUpstream(testCase)
expected=testCase.TestData.Oracle.sampler_public_members;
constructors={ ...
    "BruteForceSampler",@()radia.optuna.BruteForceSampler(Seed=7); ...
    "CmaEsSampler",@()radia.optuna.CmaEsSampler(Seed=7,NStartupTrials=99); ...
    "GPSampler",@()radia.optuna.GPSampler(Seed=7,NStartupTrials=99); ...
    "GridSampler",@()radia.optuna.GridSampler( ...
        struct("x",{{0,1}},"y",{{1,3}}),Seed=7); ...
    "NSGAIIISampler",@()radia.optuna.NSGAIIISampler( ...
        Seed=7,PopulationSize=4); ...
    "NSGAIISampler",@()radia.optuna.NSGAIISampler( ...
        Seed=7,PopulationSize=4); ...
    "PartialFixedSampler",@()radia.optuna.PartialFixedSampler( ...
        struct("x",0.25),radia.optuna.TPESampler( ...
        Seed=7,NStartupTrials=1,Multivariate=true)); ...
    "QMCSampler",@()radia.optuna.QMCSampler( ...
        Seed=7,WarnAsynchronousSeeding=false); ...
    "RandomSampler",@()radia.optuna.RandomSampler(7); ...
    "TPESampler",@()radia.optuna.TPESampler( ...
        Seed=7,NStartupTrials=1,Multivariate=true)};
for index=1:size(constructors,1)
    name=constructors{index,1};
    sampler=constructors{index,2}();
    study=publicMemberStudy();
    complete=study.get_trials("COMPLETE");
    running=study.ask();
    actual=sampler.infer_relative_search_space(study,running);
    actualKeys=sort(reshape(string({actual.name}),1,[]));
    samplerExpected=expected.samplers.(name);
    verifyEqual(testCase,actualKeys,sort(reshape(string( ...
        samplerExpected.infer_relative_search_space_keys),1,[])),name);
    verifyTrue(testCase,logical(samplerExpected.before_trial_returns_none));
    verifyWarningFree(testCase,@()sampler.before_trial(study,running),name);
    verifyTrue(testCase,logical(samplerExpected.after_trial_returns_none));
    verifyWarningFree(testCase,@()sampler.after_trial(study,complete(1), ...
        radia.optuna.TrialState.COMPLETE,complete(1).Values),name);
    directSampler=constructors{index,2}();
    directStudy=radia.optuna.Study(Sampler=directSampler,AutoSave=false);
    directTrial=directStudy.ask();
    directSpace=directSampler.infer_relative_search_space( ...
        directStudy,directTrial);
    independent=directSampler.sample_independent( ...
        directStudy,directTrial,"x", ...
        radia.optuna.FloatDistribution(0,1,Step=0.1));
    verifyEqual(testCase,independent,samplerExpected.independent_value,name);
    relative=directSampler.sample_relative( ...
        directStudy,directTrial,directSpace);
    verifyEqual(testCase,string(fieldnames(relative)), ...
        string(fieldnames(samplerExpected.relative_params)),name);
    relativeNames=string(fieldnames(relative));
    for relativeName=reshape(relativeNames,1,[])
        verifyEqual(testCase,relative.(relativeName), ...
            samplerExpected.relative_params.(relativeName),name);
    end
end

grid=radia.optuna.GridSampler(struct("x",{{0,1}}),Seed=7);
gridStudy=radia.optuna.Study(Sampler=grid,AutoSave=false);
verifyEqual(testCase,grid.is_exhausted(gridStudy), ...
    logical(expected.grid_is_exhausted_before));
for index=1:2
    trial=gridStudy.ask();
    value=trial.suggest_int("x",0,1);
    gridStudy.tell(trial,value);
end
verifyEqual(testCase,grid.is_exhausted(gridStudy), ...
    logical(expected.grid_is_exhausted_after));

verifyEqual(testCase,radia.optuna.NSGAIISampler( ...
    PopulationSize=5).population_size,double(expected.nsgaii_population_size));
verifyEqual(testCase,radia.optuna.NSGAIIISampler( ...
    PopulationSize=6).population_size,double(expected.nsgaiii_population_size));
gaMetadata=meta.class.fromName("radia.optuna.BaseGASampler");
verifyTrue(testCase,gaMetadata.Abstract);
gaConstructors={ ...
    "NSGAIISampler",@()radia.optuna.NSGAIISampler( ...
        Seed=11,PopulationSize=2); ...
    "NSGAIIISampler",@()radia.optuna.NSGAIIISampler( ...
        Seed=11,PopulationSize=2)};
for gaIndex=1:size(gaConstructors,1)
    gaName=gaConstructors{gaIndex,1};
    gaSampler=gaConstructors{gaIndex,2}();
    gaStudy=radia.optuna.Study(Sampler=gaSampler,AutoSave=false);
    for objective=[3,1,2,0]
        gaTrial=gaStudy.ask();
        gaTrial.suggest_float("x",0,1);
        gaStudy.tell(gaTrial,objective);
    end
    contract=expected.ga.(gaName);
    gaTrials=gaStudy.get_trials();
    generations=arrayfun(@(trial)gaSampler.get_trial_generation( ...
        gaStudy,trial),gaTrials);
    verifyEqual(testCase,reshape(generations,1,[]), ...
        reshape(double(contract.generations),1,[]),gaName);
    population=gaSampler.get_population(gaStudy,0);
    verifyEqual(testCase,[population.Number], ...
        reshape(double(contract.population_numbers),1,[]),gaName);
    parents=gaSampler.get_parent_population(gaStudy,1);
    verifyEqual(testCase,[parents.Number], ...
        reshape(double(contract.parent_numbers),1,[]),gaName);
    selected=gaSampler.select_parent(gaStudy,1);
    verifyEqual(testCase,[selected.Number], ...
        reshape(double(contract.selected_numbers),1,[]),gaName);
    gaSampler.population_size=3;
    verifyEqual(testCase,gaSampler.population_size, ...
        double(contract.population_size_after_set),gaName);
end

warning("off","radia:optuna:FutureWarning");
warningCleanup=onCleanup(@()warning("on","radia:optuna:FutureWarning"));
hyperopt=radia.optuna.TPESampler.hyperopt_parameters();
verifyEqual(testCase,hyperopt.consider_endpoints, ...
    logical(expected.hyperopt_parameters.consider_endpoints));
verifyEqual(testCase,hyperopt.consider_magic_clip, ...
    logical(expected.hyperopt_parameters.consider_magic_clip));
verifyEqual(testCase,hyperopt.consider_prior, ...
    logical(expected.hyperopt_parameters.consider_prior));
verifyEqual(testCase,hyperopt.n_ei_candidates, ...
    double(expected.hyperopt_parameters.n_ei_candidates));
verifyEqual(testCase,hyperopt.n_startup_trials, ...
    double(expected.hyperopt_parameters.n_startup_trials));
verifyEqual(testCase,hyperopt.prior_weight, ...
    double(expected.hyperopt_parameters.prior_weight));
verifyEqual(testCase,arrayfun(hyperopt.gamma,[0,1,16,10000]), ...
    reshape(double(expected.hyperopt_parameters.gamma),1,[]));
weightCounts=[0,3,27];
for index=1:numel(weightCounts)
    verifyEqual(testCase,hyperopt.weights(weightCounts(index)), ...
        reshape(double(expected.hyperopt_parameters.weights{index}),1,[]));
end
clear warningCleanup

crossovers={ ...
    "blxalpha",radia.optuna.nsgaii.BLXAlphaCrossover(); ...
    "sbx",radia.optuna.nsgaii.SBXCrossover(); ...
    "spx",radia.optuna.nsgaii.SPXCrossover(); ...
    "undx",radia.optuna.nsgaii.UNDXCrossover(); ...
    "uniform",radia.optuna.nsgaii.UniformCrossover(); ...
    "vsbx",radia.optuna.nsgaii.VSBXCrossover()};
for index=1:size(crossovers,1)
    verifyEqual(testCase,crossovers{index,2}.n_parents, ...
        double(expected.crossover_n_parents.(crossovers{index,1})));
end
verifyEqual(testCase,string(expected.base_crossover_instantiation_error), ...
    "TypeError");
metadata=meta.class.fromName("radia.optuna.nsgaii.BaseCrossover");
verifyTrue(testCase,metadata.Abstract);

fixed=radia.optuna.FixedTrial(struct("x",0.5));
verifyEqual(testCase,~ismissing(fixed.datetime_start()), ...
    logical(expected.fixed_trial_datetime_start_is_not_none));
relativeStudy=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=7,NStartupTrials=0,Multivariate=true),AutoSave=false);
complete=relativeStudy.ask();
complete.suggest_float("x",0,1,Step=0.1);
complete.suggest_int("y",1,5);
relativeStudy.tell(complete,0.2);
running=relativeStudy.ask();
verifyEqual(testCase,~ismissing(running.datetime_start()), ...
    logical(expected.trial_datetime_start_is_not_none));
verifyEqual(testCase,sort(reshape(string(fieldnames( ...
    running.relative_params)),1,[])),sort(reshape(string( ...
    expected.trial_relative_params_keys),1,[])));

names=sort(reshape(string(fieldnames(expected.mapped_namespaces)),1,[]));
verifyEqual(testCase,names,sort(["distributions","exceptions", ...
    "importance","pruners","samplers","search_space","storages", ...
    "study","trial"]));
for name=reshape(names,1,[])
    verifyEqual(testCase,string(expected.mapped_namespaces.(name)), ...
        "optuna."+name);
end
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

function testExceptionPublicContractMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.exceptions;
names=string(fieldnames(expected));
for name=reshape(names,1,[])
    contract=expected.(name);
    constructors={cell(1,0),{"oracle message"},{"a","b"}};
    for index=1:numel(constructors)
        arguments=constructors{index};
        exception=feval("radia.optuna."+name,arguments{:});
        item=contract.cases(index);
        verifyEqual(testCase,exception.message,string(item.message),name);
        verifyEqual(testCase,string(exception.args), ...
            reshape(string(item.args),1,[]),name);
        exception.add_note("oracle note");
        verifyEqual(testCase,exception.notes, ...
            reshape(string(item.notes),1,[]),name);
        verifyTrue(testCase,exception.with_traceback([])==exception,name);
        try
            throw(exception);
            verifyFail(testCase,"Optuna exception was not thrown: "+name);
        catch caught
            verifyEqual(testCase,string(caught.identifier), ...
                "radia:optuna:"+name,name);
            verifyEqual(testCase,string(caught.message),string(item.message),name);
        end
    end
    verifyEqual(testCase,isa(feval("radia.optuna."+name), ...
        "radia.optuna.OptunaError"),logical(contract.is_optuna_error),name);
    verifyEqual(testCase,name=="ExperimentalWarning", ...
        logical(contract.is_warning),name);
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

function testTrialsDataframeMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.trials_dataframe;
distributions=struct( ...
    "x",radia.optuna.FloatDistribution(0,1), ...
    "mode",radia.optuna.CategoricalDistribution(["A","B"]));
study=radia.optuna.Study(AutoSave=false);
study.add_trial(radia.optuna.create_trial(value=1.25, ...
    params=struct("x",0.25,"mode","A"), ...
    distributions=distributions, ...
    user_attrs=struct("owner","lab"), ...
    system_attrs=struct("origin","oracle")));
intermediate=table([0;2],[3;1], ...
    repmat(datetime("now",TimeZone="local"),2,1), ...
    VariableNames=["Step","Value","Timestamp"]);
study.add_trial(radia.optuna.create_trial(state="PRUNED", ...
    params=struct("x",0.75,"mode","B"), ...
    distributions=distributions,intermediate_values=intermediate, ...
    user_attrs=struct("owner","mdx"), ...
    system_attrs=struct("origin","oracle")));
study.add_trial(radia.optuna.create_trial(state="FAIL", ...
    params=struct("x",0.5,"mode","A"), ...
    distributions=distributions, ...
    user_attrs=struct("owner","lab"), ...
    system_attrs=struct("origin","imported")));

attrs=reshape(string(expected.attrs),1,[]);
flat=study.trials_dataframe(attrs=attrs);
verifyEqual(testCase,string(flat.Properties.VariableNames), ...
    reshape(string(expected.single.flat_columns),1,[]));
verifyEqual(testCase, ...
    reshape(string(flat.Properties.UserData.flat_columns),1,[]), ...
    reshape(string(expected.single.flat_columns),1,[]));
verifyEqual(testCase,flat.number, ...
    reshape(double(expected.single.values.number),[],1));
verifyEqual(testCase,flat.value, ...
    reshape(double(expected.single.values.value),[],1));
verifyEqual(testCase,string(flat.params_mode), ...
    reshape(string(expected.single.values.params_mode),[],1));
verifyEqual(testCase,flat.params_x, ...
    reshape(double(expected.single.values.params_x),[],1));
verifyEqual(testCase,string(flat.user_attrs_owner), ...
    reshape(string(expected.single.values.user_attrs_owner),[],1));
verifyEqual(testCase,string(flat.system_attrs_origin), ...
    reshape(string(expected.single.values.system_attrs_origin),[],1));
verifyEqual(testCase,string(flat.state), ...
    reshape(string(expected.single.values.state),[],1));

multi=study.trials_dataframe(attrs=attrs,multi_index=true);
levels=multi.Properties.UserData.column_levels;
verifyTrue(testCase,multi.Properties.UserData.multi_index);
verifyEqual(testCase,levels(:,1), ...
    reshape(string({expected.single.multi_columns.top}),[],1));
verifyEqual(testCase,levels(:,2), ...
    reshape(string({expected.single.multi_columns.sub}),[],1));

metricStudy=radia.optuna.Study( ...
    Directions=["minimize","maximize"],AutoSave=false);
metricStudy.set_metric_names(["loss","gain"]);
metricStudy.add_trial(radia.optuna.create_trial(values=[1.5,2.5], ...
    params=struct("x",0.5), ...
    distributions=struct( ...
        "x",radia.optuna.FloatDistribution(0,1))));
metricFrame=metricStudy.trials_dataframe( ...
    attrs=["number","value","params","state"],multi_index=true);
verifyEqual(testCase,string(metricFrame.Properties.VariableNames), ...
    reshape(string(expected.metric_names.flat_columns),1,[]));
verifyEqual(testCase,metricFrame.values_gain, ...
    double(expected.metric_names.values.values_gain));
verifyEqual(testCase,metricFrame.values_loss, ...
    double(expected.metric_names.values.values_loss));
metricLevels=metricFrame.Properties.UserData.column_levels;
verifyEqual(testCase,metricLevels(:,1), ...
    reshape(string({expected.metric_names.multi_columns.top}),[],1));
verifyEqual(testCase,metricLevels(:,2), ...
    reshape(string({expected.metric_names.multi_columns.sub}),[],1));

singleMetricStudy=radia.optuna.Study(AutoSave=false);
singleMetricStudy.set_metric_names("loss");
singleMetricStudy.add_trial(radia.optuna.create_trial(value=1.2));
singleMetricFrame=singleMetricStudy.trials_dataframe( ...
    attrs=["number","value","state"],multi_index=true);
verifyEqual(testCase,string(singleMetricFrame.Properties.VariableNames), ...
    reshape(string(expected.single_metric_name.flat_columns),1,[]));
verifyEqual(testCase,singleMetricFrame.value_loss, ...
    double(expected.single_metric_name.values.value_loss));
singleMetricLevels=singleMetricFrame.Properties.UserData.column_levels;
verifyEqual(testCase,singleMetricLevels(:,1), ...
    reshape(string({expected.single_metric_name.multi_columns.top}),[],1));
verifyEqual(testCase,singleMetricLevels(:,2), ...
    reshape(string({expected.single_metric_name.multi_columns.sub}),[],1));

emptyStudy=radia.optuna.Study(AutoSave=false);
emptyFrame=emptyStudy.trials_dataframe();
verifyEqual(testCase,width(emptyFrame),double(expected.empty_column_count));
verifyFalse(testCase,logical(expected.multi_index_default));
verifyEqual(testCase,string(expected.errors.unknown_attr),"AttributeError");
verifyError(testCase,@()study.trials_dataframe(attrs="not_an_attr"), ...
    "radia:optuna:TrialsDataframeAttribute");
verifyEqual(testCase,string(expected.errors.empty_attrs),"TypeError");
verifyError(testCase,@()study.trials_dataframe(attrs=strings(1,0)), ...
    "radia:optuna:TrialsDataframeAttrs");
end

function testParameterImportancesMatchUpstream(testCase)
expected=testCase.TestData.Oracle.importance;
study=importanceStudy(expected);
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

function testImportanceEvaluatorPublicMembersMatchUpstream(testCase)
expected=testCase.TestData.Oracle.importance;
contract=expected.public_evaluators;
verifyEqual(testCase,string(contract.base_construction_error),"TypeError");
verifyError(testCase,@()radia.optuna.BaseImportanceEvaluator(), ...
    "MATLAB:class:abstract");
study=importanceStudy(expected);
evaluators={ ...
    "fanova",radia.optuna.FanovaImportanceEvaluator( ...
        n_trees=16,max_depth=16,seed=97); ...
    "mdi",radia.optuna.MeanDecreaseImpurityImportanceEvaluator( ...
        n_trees=16,max_depth=16,seed=97); ...
    "ped_anova",radia.optuna.PedAnovaImportanceEvaluator( ...
        target_quantile=0.25,region_quantile=1)};
for index=1:size(evaluators,1)
    name=evaluators{index,1};
    actual=evaluators{index,2}.evaluate(study);
    verifyImportanceTable(testCase,actual,contract.(name));
end
target=radia.optuna.FanovaImportanceEvaluator( ...
    n_trees=16,max_depth=16,seed=97).evaluate( ...
    study,target=@(trial)trial.Params.y);
verifyImportanceTable(testCase,target,contract.fanova_target);
end

function testTerminationContractsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.terminator;
for name=["BaseErrorEvaluator","BaseImprovementEvaluator","BaseTerminator"]
    metadata=meta.class.fromName("radia.optuna."+name);
    verifyTrue(testCase,metadata.Abstract,name);
end

cvStudy=radia.optuna.Study(Directions="maximize",AutoSave=false);
rows=double(expected.cross_validation_rows);
for index=1:size(rows,1)
    trial=cvStudy.ask();
    radia.optuna.report_cross_validation_scores(trial,rows(index,:));
    cvStudy.tell(trial,mean(rows(index,:)));
end
cvEvaluator=radia.optuna.CrossValidationErrorEvaluator();
verifyEqual(testCase,cvEvaluator.evaluate( ...
    cvStudy.get_trials(),cvStudy.direction()), ...
    double(expected.cross_validation_error),AbsTol=1e-15);
trial=cvStudy.ask();
verifyError(testCase,@()radia.optuna.report_cross_validation_scores( ...
    trial,1),"radia:optuna:CrossValidationScores");
cvStudy.tell(trial,State="FAIL");

staticEvaluator=radia.optuna.StaticErrorEvaluator(1.25);
verifyEqual(testCase,staticEvaluator.evaluate([], ...
    radia.optuna.StudyDirection.MINIMIZE),double(expected.static_error));

medianEvaluator=radia.optuna.MedianErrorEvaluator( ...
    radia.optuna.BestValueStagnationEvaluator(5), ...
    WarmUpTrials=1,NInitialTrials=3,ThresholdRatio=0.1);
medianValues=[5,4,4.5,4.2];
medianTrials=arrayfun(@(index)radia.optuna.create_trial( ...
    value=medianValues(index)),1:4);
verifyEqual(testCase,medianEvaluator.evaluate(medianTrials,"minimize"), ...
    double(expected.median_error));
verifyEqual(testCase,medianEvaluator.evaluate( ...
    radia.optuna.create_trial(value=100),"minimize"), ...
    double(expected.median_cached));

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

function testArtifactPublicContractMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.artifacts;
tempDirectory=string(tempname("C:\temp"));
mkdir(tempDirectory);
cleanup=onCleanup(@()rmdir(tempDirectory,"s"));
source=fullfile(tempDirectory,"oracle.txt");
sourceBytes=uint8([97,114,116,105,102,97,99,116,45,98,121,116,101,115,0,255]);
radia.optuna.internal.ArtifactIO.writeFile(source,sourceBytes,false);
store=radia.optuna.FileSystemArtifactStore(tempDirectory);
study=radia.optuna.Study(AutoSave=false);
artifactId=radia.optuna.upload_artifact( ...
    artifact_store=store,file_path=source,study_or_trial=study);
verifyEqual(testCase,strlength(artifactId),strlength(string(expected.artifact_id)));
metadata=radia.optuna.get_all_artifact_meta(study);
verifyEqual(testCase,numel(metadata),1);
verifyEqual(testCase,metadata.filename,string(expected.metadata.filename));
verifyEqual(testCase,metadata.mimetype,string(expected.metadata.mimetype));
verifyTrue(testCase,ismissing(metadata.encoding));
destination=fullfile(tempDirectory,"downloaded.txt");
radia.optuna.download_artifact(artifact_store=store, ...
    file_path=destination,artifact_id=artifactId);
verifyEqual(testCase,radia.optuna.internal.ArtifactIO.readFile(destination), ...
    reshape(uint8(expected.downloaded),[],1));
verifyEqual(testCase,string(expected.existing_download_error),"FileExistsError");
verifyError(testCase,@()radia.optuna.download_artifact( ...
    artifact_store=store,file_path=destination,artifact_id=artifactId), ...
    "radia:optuna:ArtifactFileExists");
verifyEqual(testCase,string(expected.traversal_error),"ValueError");
verifyError(testCase,@()store.open_reader("../outside"), ...
    "radia:optuna:ArtifactId");

backoff=radia.optuna.Backoff(store,MaxRetries=2, ...
    MinDelay=1e-9,MaxDelay=2e-9);
backoff.write("backoff",uint8(expected.backoff_body));
verifyEqual(testCase,backoff.open_reader("backoff"), ...
    reshape(uint8(expected.backoff_body),[],1));
verifyEqual(testCase,string(expected.backoff_remove_error),"ArtifactNotFound");
verifyError(testCase,@()backoff.remove("backoff"), ...
    "radia:optuna:ArtifactNotFound");

botoClient=struct( ...
    "get_object",@(~,~)reshape(uint8(expected.boto_open),[],1), ...
    "upload_fileobj",@(~,~,~)[],"delete_object",@(~,~)[]);
boto=radia.optuna.Boto3ArtifactStore("bucket",botoClient);
verifyEqual(testCase,boto.open_reader("cloud"), ...
    reshape(uint8(expected.boto_open),[],1));
verifyWarningFree(testCase,@()boto.write("cloud",uint8(expected.boto_written)));
verifyWarningFree(testCase,@()boto.remove("cloud"));

gcsClient=struct( ...
    "get_blob",@(~,~)reshape(uint8(expected.gcs_open),[],1), ...
    "upload_blob",@(~,~,~)[],"delete_blob",@(~,~)[]);
gcs=radia.optuna.GCSArtifactStore("bucket",gcsClient);
verifyEqual(testCase,gcs.open_reader("cloud"), ...
    reshape(uint8(expected.gcs_open),[],1));
verifyWarningFree(testCase,@()gcs.write("cloud",uint8(expected.gcs_written)));
verifyWarningFree(testCase,@()gcs.remove("cloud"));
clear cleanup
end

function testInMemoryStorageMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.storage;
metadata=meta.class.fromName("radia.optuna.BaseStorage");
verifyEqual(testCase,metadata.Abstract,logical(expected.base_is_abstract));
storage=radia.optuna.InMemoryStorage();
studyId=storage.create_new_study( ...
    radia.optuna.StudyDirection.MINIMIZE,"memory-oracle");
verifyEqual(testCase,studyId,double(expected.study_id));
storage.set_study_user_attr(studyId,"owner","oracle");
storage.set_study_system_attr(studyId,"revision",4);
trialId=storage.create_new_trial(studyId);
verifyEqual(testCase,trialId,double(expected.trial_id));
storage.set_trial_param(trialId,"x",0.5, ...
    radia.optuna.FloatDistribution(0,1));
storage.set_trial_user_attr(trialId,"label","first");
storage.set_trial_system_attr(trialId,"worker",7);
storage.set_trial_intermediate_value(trialId,2,3.5);
verifyEqual(testCase,storage.get_trial(trialId).State, ...
    string(expected.running_state));
verifyEqual(testCase,storage.set_trial_state_values( ...
    trialId,radia.optuna.TrialState.COMPLETE,1.25), ...
    logical(expected.completed));
frozen=storage.get_trial(trialId);
verifyEqual(testCase,frozen.State,string(expected.trial_state));
verifyEqual(testCase,frozen.Value,double(expected.trial_value));
verifyEqual(testCase,storage.get_trial_param(trialId,"x"), ...
    double(expected.param_internal));
verifyEqual(testCase,storage.get_trial_params(trialId),expected.params);
trialUserAttrs=storage.get_trial_user_attrs(trialId);
verifyEqual(testCase,string(trialUserAttrs.label), ...
    string(expected.trial_user_attrs.label));
verifyEqual(testCase,storage.get_trial_system_attrs(trialId), ...
    expected.trial_system_attrs);
verifyEqual(testCase,string(expected.finished_update_error), ...
    "UpdateFinishedTrialError");
verifyError(testCase,@()storage.set_trial_user_attr( ...
    trialId,"late",true),"radia:optuna:UpdateFinishedTrialError");

failedId=storage.create_new_trial(studyId);
verifyEqual(testCase,failedId,double(expected.failed_id));
storage.set_trial_state_values(failedId,radia.optuna.TrialState.FAIL);
templateId=storage.create_new_trial(studyId, ...
    radia.optuna.create_trial(value=0.25));
verifyEqual(testCase,templateId,double(expected.template_id));
verifyEqual(testCase,storage.get_best_trial(studyId).Number, ...
    double(expected.best_number));
verifyEqual(testCase,storage.get_n_trials(studyId),double(expected.n_trials));
verifyEqual(testCase,storage.get_n_trials( ...
    studyId,radia.optuna.TrialState.COMPLETE),double(expected.n_complete));
complete=storage.get_all_trials( ...
    studyId,true,radia.optuna.TrialState.COMPLETE);
verifyEqual(testCase,[complete.Number], ...
    reshape(double(expected.complete_numbers),1,[]));
verifyEqual(testCase,storage.get_trial_id_from_study_id_trial_number( ...
    studyId,frozen.Number),double(expected.trial_lookup_id));
verifyEqual(testCase,storage.get_trial_number_from_id(trialId), ...
    double(expected.trial_number));
verifyEqual(testCase,string(storage.get_study_directions(studyId)), ...
    reshape(string(expected.directions),1,[]));
verifyEqual(testCase,storage.get_study_name_from_id(studyId), ...
    string(expected.study_name));
verifyEqual(testCase,storage.get_study_id_from_name("memory-oracle"),studyId);
studyUserAttrs=storage.get_study_user_attrs(studyId);
verifyEqual(testCase,string(studyUserAttrs.owner), ...
    string(expected.study_user_attrs.owner));
verifyEqual(testCase,storage.get_study_system_attrs(studyId), ...
    expected.study_system_attrs);
summaries=storage.get_all_studies();
verifyEqual(testCase,numel(summaries),1);
verifyEqual(testCase,summaries(1).study_name,string(expected.summary.name));
verifyEqual(testCase,string(expected.duplicate_error),"DuplicatedStudyError");
verifyError(testCase,@()storage.create_new_study( ...
    radia.optuna.StudyDirection.MINIMIZE,"memory-oracle"), ...
    "radia:optuna:DuplicatedStudyError");
verifyTrue(testCase,logical(expected.remove_session_is_none));
verifyWarningFree(testCase,@()storage.remove_session());
storage.delete_study(studyId);
verifyError(testCase,@()storage.get_study_name_from_id(studyId), ...
    "radia:optuna:StudyNotFound");
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

function testDistributionPublicMembersMatchUpstream(testCase)
expected=testCase.TestData.Oracle.distribution_public_members;
verifyEqual(testCase,string(expected.base_construction_error),"TypeError");
verifyError(testCase,@()radia.optuna.BaseDistribution(), ...
    "MATLAB:class:abstract");

warning("off","radia:optuna:FutureWarning");
cleanup=onCleanup(@()warning("on","radia:optuna:FutureWarning"));
instances=struct( ...
    "FloatDistribution",radia.optuna.FloatDistribution(1,1), ...
    "IntDistribution",radia.optuna.IntDistribution(1,1), ...
    "CategoricalDistribution", ...
        radia.optuna.CategoricalDistribution({"A","B",3}), ...
    "UniformDistribution",radia.optuna.UniformDistribution(1,1), ...
    "LogUniformDistribution",radia.optuna.LogUniformDistribution(1,2), ...
    "DiscreteUniformDistribution", ...
        radia.optuna.DiscreteUniformDistribution(0,1,0.2), ...
    "IntUniformDistribution",radia.optuna.IntUniformDistribution(1,3,1), ...
    "IntLogUniformDistribution", ...
        radia.optuna.IntLogUniformDistribution(1,3,1));
for name=reshape(string(fieldnames(instances)),1,[])
    verifyEqual(testCase,instances.(name).single(), ...
        logical(expected.single.(name)));
end

categorical=instances.CategoricalDistribution;
verifyEqual(testCase,categorical.to_internal_repr("B"), ...
    double(expected.categorical_internal));
verifyEqual(testCase,categorical.to_external_repr(2), ...
    expected.categorical_external);
verifyEqual(testCase,instances.FloatDistribution.to_internal_repr("1.25"), ...
    double(expected.float_internal));
verifyEqual(testCase,instances.FloatDistribution.to_external_repr(1.25), ...
    double(expected.float_external));
verifyEqual(testCase,instances.IntDistribution.to_internal_repr("3"), ...
    double(expected.int_internal));
verifyEqual(testCase,instances.IntDistribution.to_external_repr(3.9), ...
    double(expected.int_external));
verifyEqual(testCase,instances.DiscreteUniformDistribution.q, ...
    double(expected.discrete_q));
verifyEqual(testCase,radia.optuna.DISTRIBUTION_CLASSES(), ...
    reshape(string(expected.distribution_classes),1,[]));
verifyNotEmpty(testCase,string(expected.categorical_choice_type));
verifyNotEmpty(testCase,radia.optuna.CategoricalChoiceType());
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

function testGroupDecomposedSearchSpaceMatchesUpstream(testCase)
expected=testCase.TestData.Oracle.search_space.group;
base=struct( ...
    "x",radia.optuna.FloatDistribution(0,1), ...
    "fixed",radia.optuna.FloatDistribution(2,2));
group=radia.optuna.SearchSpaceGroup();
group.add_distributions(struct("x",base.x,"y",base.fixed));
group.add_distributions(struct( ...
    "x",base.x,"z",radia.optuna.IntDistribution(1,3)));
verifyEqual(testCase,groupSignatures(group.search_spaces()), ...
    reshape(string(expected.direct_signatures),[],1));

study=radia.optuna.Study(AutoSave=false);
study.add_trials([ ...
    radia.optuna.create_trial(value=1,params=struct("x",0.5,"fixed",2), ...
        distributions=struct("x",base.x,"fixed",base.fixed)); ...
    radia.optuna.create_trial(value=1,params=struct("x",0.5,"z",1), ...
        distributions=struct("x",base.x,"z", ...
        radia.optuna.IntDistribution(1,3)))]);
calculator=radia.optuna.GroupDecomposedSearchSpace();
calculated=calculator.calculate(study);
verifyEqual(testCase,groupSignatures(calculated.search_spaces()), ...
    reshape(string(expected.calculated_signatures),[],1));
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

enumCases={ ...
    radia.optuna.StudyDirection.MAXIMIZE,expected.integer_api.study_direction; ...
    radia.optuna.TrialState.WAITING,expected.integer_api.trial_state};
for index=1:size(enumCases,1)
    value=enumCases{index,1};
    contract=enumCases{index,2};
    verifyEqual(testCase,value.as_integer_ratio(), ...
        reshape(double(contract.as_integer_ratio),1,[]));
    verifyEqual(testCase,value.bit_count(),double(contract.bit_count));
    verifyEqual(testCase,value.bit_length(),double(contract.bit_length));
    verifyEqual(testCase,value.conjugate(),double(contract.conjugate));
    verifyEqual(testCase,value.denominator(),double(contract.denominator));
    verifyEqual(testCase,value.imag(),double(contract.imag));
    verifyEqual(testCase,value.is_integer(),logical(contract.is_integer));
    verifyEqual(testCase,value.name(),string(contract.name));
    verifyEqual(testCase,value.numerator(),double(contract.numerator));
    verifyEqual(testCase,value.real(),double(contract.real));
    verifyEqual(testCase,value.to_bytes(2,"big"), ...
        reshape(uint8(contract.to_bytes),1,[]));
    verifyEqual(testCase,value.value(),double(contract.value));
end
verifyEqual(testCase,radia.optuna.StudyDirection.from_bytes( ...
    uint8([0,2]),"big").name(), ...
    string(expected.integer_api.study_direction.from_bytes_name));
verifyEqual(testCase,radia.optuna.TrialState.from_bytes( ...
    uint8([0,4]),"big").name(), ...
    string(expected.integer_api.trial_state.from_bytes_name));

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
verifyEqual(testCase,percentile.Pruner.prune(percentile,trial), ...
    logical(expected.percentile_minimize));

maximize=radia.optuna.Study(Directions="maximize",Pruner= ...
    radia.optuna.PercentilePruner(50,NStartupTrials=0),AutoSave=false);
addCompletedTrial(maximize,0,1); addCompletedTrial(maximize,0,3);
trial=maximize.ask(); trial.report(1.5,0);
verifyEqual(testCase,maximize.Pruner.prune(maximize,trial), ...
    logical(expected.percentile_maximize));

median=radia.optuna.Study(Pruner=radia.optuna.MedianPruner( ...
    NStartupTrials=0),AutoSave=false);
addCompletedTrial(median,0,1); addCompletedTrial(median,0,3);
medianTrial=median.ask(); medianTrial.report(5,0);
verifyEqual(testCase,median.Pruner.prune(median,medianTrial), ...
    logical(expected.median));

threshold=radia.optuna.Study(Pruner=radia.optuna.ThresholdPruner( ...
    Lower=0,Upper=10),AutoSave=false);
trial=threshold.ask(); trial.report(5,0);
actual(1)=threshold.Pruner.prune(threshold,trial);
trial.report(11,1); actual(2)=threshold.Pruner.prune(threshold,trial);
nanTrial=threshold.ask(); nanTrial.report(NaN,0);
actual(3)=threshold.Pruner.prune(threshold,nanTrial);
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
verifyEqual(testCase,[patient.Pruner.prune(patient,trial), ...
    wrapped.Pruner.prune(wrapped,wrappedTrial)], ...
    reshape(logical(expected.patient),1,[]));

halvingPruner=radia.optuna.SuccessiveHalvingPruner( ...
    MinResource=1,ReductionFactor=2);
halving=radia.optuna.Study(Pruner=halvingPruner,AutoSave=false);
first=halving.ask(); first.report(1,1);
firstDecision=halvingPruner.prune(halving,first); halving.tell(first,1);
second=halving.ask(); second.report(2,1);
secondDecision=halvingPruner.prune(halving,second);
verifyEqual(testCase,[firstDecision,secondDecision], ...
    reshape(logical(expected.successive_halving.decisions),1,[]));
verifyEqual(testCase,[first.SystemAttrs.completed_rung_0, ...
    second.SystemAttrs.completed_rung_0], ...
    reshape(double(expected.successive_halving.rung_values),1,[]));
bootstrap=radia.optuna.Study(Pruner= ...
    radia.optuna.SuccessiveHalvingPruner(MinResource=1, ...
    ReductionFactor=2,BootstrapCount=1),AutoSave=false);
bootstrapTrial=bootstrap.ask(); bootstrapTrial.report(1,1);
verifyEqual(testCase,bootstrap.Pruner.prune(bootstrap,bootstrapTrial), ...
    logical(expected.successive_halving.bootstrap));

hyperbandPruner=radia.optuna.HyperbandPruner(MinResource=1, ...
    MaxResource=9,ReductionFactor=3);
hyperband=radia.optuna.Study(Name="hb",Pruner=hyperbandPruner,AutoSave=false);
hyperbandTrial=hyperband.ask(); hyperbandTrial.report(1,0);
verifyEqual(testCase,hyperbandPruner.prune(hyperband,hyperbandTrial), ...
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
verifyEqual(testCase,[wilcoxon.Pruner.prune(wilcoxon,wilcoxonTrial), ...
    wilcoxon.Pruner.prune(wilcoxon,nonfinite)], ...
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

function study=publicMemberStudy()
study=radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(91),AutoSave=false);
for value=[0.2,0.1]
    trial=study.ask();
    trial.suggest_float("x",0,1,Step=0.1);
    trial.suggest_int("y",1,5);
    study.tell(trial,value);
end
end

function value=failingConstraintCallback()
value=[]; %#ok<NASGU>
error("radia:test:ConstraintCallback","constraint callback failed");
end

function value=readSummaryDirection(summary)
value=summary.direction;
end

function study=importanceStudy(expected)
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
end

function verifyImportanceTable(testCase,actual,expected)
expectedNames=reshape(string(expected.parameter_order),[],1);
expectedValues=reshape(double(expected.values),[],1);
verifyEqual(testCase,sort(actual.Parameter),sort(expectedNames));
for index=1:numel(expectedNames)
    row=find(actual.Parameter==expectedNames(index),1);
    verifyNotEmpty(testCase,row);
    verifyEqual(testCase,actual.Importance(row),expectedValues(index),AbsTol=0);
end
end

function signatures=groupSignatures(spaces)
signatures=strings(numel(spaces),1);
for index=1:numel(spaces)
    signatures(index)=strjoin(sort(string(spaces{index}.keys)),",");
end
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
