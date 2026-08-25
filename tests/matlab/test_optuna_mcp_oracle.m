function tests=test_optuna_mcp_oracle
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
testCase.TestData.RemovePath=~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemovePath, addpath(matlabDirectory); end
testCase.TestData.MatlabDirectory=matlabDirectory;
fixture=fullfile(root,"tests","matlab","fixtures", ...
    "optuna49_mcp_oracle.json");
testCase.TestData.Oracle=jsondecode(fileread(fixture));
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testOfficialMCPProvenanceAndToolContract(testCase)
oracle=testCase.TestData.Oracle;
verifyEqual(testCase,string(oracle.schema), ...
    "radia.test.optuna-upstream-mcp-oracle.v1");
verifyEqual(testCase,string(oracle.transport),"stdio");
verifyEqual(testCase,string(oracle.optuna_version),"4.9.0");
verifyEqual(testCase,string(oracle.optuna_mcp_version),"0.2.0");
verifyEqual(testCase,string(oracle.mcp_server_name),"Optuna");
verifyNotEmpty(testCase,string(oracle.mcp_server_reported_version));

tools=string(oracle.tools);
required=["create_study","set_sampler","ask","tell", ...
    "set_trial_user_attr","get_trial_user_attrs","get_trials", ...
    "best_trial","best_trials","add_trial","add_trials", ...
    "set_metric_names","get_metric_names","get_directions"];
verifyTrue(testCase,all(ismember(required,tools)));
verifyFalse(testCase,logical(oracle.sampler_seed_supported));
verifyEqual(testCase,string(oracle.set_sampler_arguments),"name");
end

function testSingleObjectiveStudyTrialContract(testCase)
expected=testCase.TestData.Oracle.single;
study=radia.optuna.Study(Name=string(expected.study_name), ...
    Directions=string(expected.directions.directions), ...
    Sampler=radia.optuna.RandomSampler(0),AutoSave=false);
verifyEqual(testCase,study.Name,string(expected.create_study.study_name));
verifyEqual(testCase,string(class(study.Sampler)), ...
    "radia.optuna."+string(expected.set_sampler.sampler_name));

fixed=struct( ...
    "x",radia.optuna.FloatDistribution(1.25,1.25), ...
    "n",radia.optuna.IntDistribution(3,3), ...
    "mode",radia.optuna.CategoricalDistribution("A"));
trial=study.ask(fixed);
verifyEqual(testCase,trial.Number,double(expected.ask.trial_number));
verifyEqual(testCase,trial.Params.x,double(expected.ask.params.x));
verifyEqual(testCase,trial.Params.n,double(expected.ask.params.n));
verifyEqual(testCase,string(trial.Params.mode),string(expected.ask.params.mode));

trial.set_user_attr("source", ...
    string(expected.user_attrs.user_attrs.source));
verifyEqual(testCase,string(trial.user_attrs().source), ...
    string(expected.user_attrs.user_attrs.source));
reported=study.tell(trial,double(expected.tell.values));
verifyEqual(testCase,reported.Number,double(expected.tell.trial_number));
verifyEqual(testCase,reported.Values, ...
    reshape(double(expected.tell.values),1,[]));

distribution=struct("x",radia.optuna.FloatDistribution(-1,2));
study.add_trial(radia.optuna.createTrial(State="COMPLETE",Values=0.5, ...
    Params=struct("x",-0.5),Distributions=distribution, ...
    UserAttrs=struct("source","archive"), ...
    SystemAttrs=struct("origin","official-mcp")));
imports(1)=radia.optuna.createTrial(State="PRUNED", ...
    Params=struct("x",0),Distributions=distribution);
imports(2)=radia.optuna.createTrial(State="FAIL", ...
    Params=struct("x",0.75),Distributions=distribution);
study.add_trials(imports);

summary=expected.trials;
verifyEqual(testCase,study.TrialTable.TrialNumber, ...
    reshape(double(summary.numbers),[],1));
verifyEqual(testCase,study.TrialTable.State, ...
    reshape(string(summary.states),[],1));
present=~isnan(study.TrialTable.Value);
verifyEqual(testCase,present,reshape(logical(summary.value_present),[],1));
verifyEqual(testCase,study.TrialTable.Value(present), ...
    reshape(double(summary.values(logical(summary.value_present))),[],1));

actualX=zeros(height(study.TrialTable),1);
for index=1:height(study.TrialTable)
    actualX(index)=study.TrialTable.Params{index}.x;
end
verifyEqual(testCase,actualX,reshape(double(summary.params_x),[],1));
verifyEqual(testCase,study.freezeTrial(0).Params.n, ...
    double(summary.params_n(1)));
verifyEqual(testCase,string(study.freezeTrial(0).Params.mode), ...
    string(summary.params_mode(1)));
verifyEqual(testCase,string(study.freezeTrial(0).UserAttrs.source), ...
    string(summary.user_attrs_source(1)));
verifyEqual(testCase,string(study.freezeTrial(1).UserAttrs.source), ...
    string(summary.user_attrs_source(2)));

best=study.best_trial();
verifyEqual(testCase,best.Number,double(expected.best_trial.trial_number));
verifyEqual(testCase,best.Values, ...
    reshape(double(expected.best_trial.values),1,[]));
verifyEqual(testCase,best.Params.x,double(expected.best_trial.params.x));
verifyEqual(testCase,string(best.UserAttrs.source), ...
    string(expected.best_trial.user_attrs.source));
verifyEqual(testCase,string(best.SystemAttrs.origin), ...
    string(expected.best_trial.system_attrs.origin));
end

function testMultiObjectiveMetricAndParetoContract(testCase)
expected=testCase.TestData.Oracle.multi;
study=radia.optuna.Study(Name=string(expected.study_name), ...
    Directions=string(expected.directions.directions),AutoSave=false);
study.set_metric_names(string(expected.set_metric_names.metric_names));
verifyEqual(testCase,study.metric_names(), ...
    reshape(string(expected.get_metric_names.metric_names),1,[]));

trial=study.ask(struct( ...
    "x",radia.optuna.FloatDistribution(0,0)));
verifyEqual(testCase,trial.Number,double(expected.ask.trial_number));
verifyEqual(testCase,trial.Params.x,double(expected.ask.params.x));
reported=study.tell(trial,reshape(double(expected.tell.values),1,[]));
verifyEqual(testCase,reported.Values, ...
    reshape(double(expected.tell.values),1,[]));

distribution=struct("x",radia.optuna.FloatDistribution(-1,2));
imports(1)=radia.optuna.createTrial(State="COMPLETE",Values=[0.5,0.5], ...
    Params=struct("x",0.5),Distributions=distribution);
imports(2)=radia.optuna.createTrial(State="COMPLETE",Values=[2,2], ...
    Params=struct("x",2),Distributions=distribution);
imports(3)=radia.optuna.createTrial(State="COMPLETE",Values=[1.5,0], ...
    Params=struct("x",1.5),Distributions=distribution);
study.add_trials(imports);

mcpFront=expected.best_trials.result;
expectedNumbers=reshape(double([mcpFront.trial_number]),[],1);
front=study.best_trials();
verifyEqual(testCase,reshape([front.Number],[],1),expectedNumbers);
for index=1:numel(front)
    verifyEqual(testCase,front(index).Values, ...
        reshape(double(mcpFront(index).values),1,[]));
    verifyEqual(testCase,front(index).Params.x, ...
        double(mcpFront(index).params.x));
end
verifyTrue(testCase,logical(expected.best_trial_is_error));
verifyError(testCase,@()study.best_trial(), ...
    "radia:optuna:MultiObjectiveBest");
end
