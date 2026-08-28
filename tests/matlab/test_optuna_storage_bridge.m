function tests=test_optuna_storage_bridge
tests=functiontests(localfunctions);
end

function testMatlabHandoffPreservesNamesAttributesAndConstraints(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
study=radia.optuna.Study(Name="bridge-demo", ...
    Sampler=radia.optuna.RandomSampler(31),AutoSave=false);
study.set_metric_names("loss");
study.set_user_attr("owner","radia");
study.set_system_attr("source_id",7);
trial=study.ask();
value=trial.suggest_float("x-1",-1,1);
trial.set_user_attr("trial-tag","t1");
trial.set_system_attr("worker-id",3);
trial.report(1.5,0);
study.recordConstraints(trial,[-0.5 0]);
study.tell(trial,value^2);

path=string(tempname)+".json";
cleanup=onCleanup(@()deleteIfPresent(path));
payload=radia.optuna.export_study(study,Path=path);
verifyEqual(testCase,string(payload.schema), ...
    "radia.optuna.study-export.v1");
verifyTrue(testCase,isfile(path));
restored=radia.optuna.import_study(payload, ...
    Sampler=radia.optuna.RandomSampler(31),AutoSave=false);

verifyEqual(testCase,restored.ParamTable.Name,"x-1");
verifyEqual(testCase,restored.UserAttrTable.Name,"trial-tag");
verifyEqual(testCase,restored.SystemAttrTable.Name,"worker-id");
verifyEqual(testCase,restored.constraintsForTrial(0),[-0.5 0]);
verifyEqual(testCase,restored.MetricNames,"loss");
verifyEqual(testCase,string(restored.UserAttrs.owner),"radia");
verifyEqual(testCase,restored.SystemAttrs.source_id,7);
clear cleanup
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
