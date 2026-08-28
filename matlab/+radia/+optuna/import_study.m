function study=import_study(payload,options)
%IMPORT_STUDY Rebuild a study from radia.optuna.study-export.v1.

arguments
    payload
    options.Name (1,1) string = ""
    options.StoragePath (1,1) string = ""
    options.Sampler = []
    options.AutoSave (1,1) logical = true
end

if isstring(payload) || ischar(payload)
    path=string(payload);
    if ~isfile(path)
        error("radia:optuna:ImportStudy", ...
            "Export file '%s' does not exist.",path);
    end
    payload=jsondecode(fileread(path));
end
if ~isstruct(payload) || ~isscalar(payload) || ...
        ~isfield(payload,"schema") || ...
        string(payload.schema)~="radia.optuna.study-export.v1"
    error("radia:optuna:ImportStudy", ...
        "Expected one radia.optuna.study-export.v1 payload.");
end

name=options.Name;
if strlength(name)==0 && isfield(payload,"study_name")
    name=string(payload.study_name);
end
directions=reshape(string(payload.directions),1,[]);
study=radia.optuna.Study(Name=name,Directions=directions, ...
    Sampler=options.Sampler,StoragePath=options.StoragePath, ...
    AutoSave=options.AutoSave);
if isfield(payload,"metric_names") && ~isempty(payload.metric_names)
    study.set_metric_names(reshape(string(payload.metric_names),1,[]));
end
applyStudyAttributes(study,payload,"user_attrs",@study.set_user_attr);
applyStudyAttributes(study,payload,"system_attrs",@study.set_system_attr);

trials=normalizeList(payload,"trials");
numbers=cellfun(@(trial)double(trial.number),trials);
[~,order]=sort(numbers);
for index=reshape(order,1,[])
    [frozen,originalNames,userAttrNames,systemAttrNames]=frozenFromRecord( ...
        trials{index},numel(directions));
    study.addTrial(frozen,OriginalNames=originalNames, ...
        OriginalUserAttrNames=userAttrNames, ...
        OriginalSystemAttrNames=systemAttrNames);
end
end

function [frozen,originalNames,userAttrNames,systemAttrNames]= ...
        frozenFromRecord(record,objectiveCount)
entries=normalizeList(record,"params");
params=struct();
distributions=struct();
originalNames=struct();
claimed=radia.optuna.Trial.claimKeys( ...
    cellfun(@(item)string(item.name),entries));
for index=1:numel(entries)
    item=entries{index};
    key=char(claimed(index));
    params.(key)=item.value;
    distributions.(key)=radia.optuna.json_to_distribution( ...
        string(item.distribution));
    originalNames.(key)=char(string(item.name));
end

state=string(record.state);
values=NaN;
if isfield(record,"values") && ~isempty(record.values)
    values=reshape(double(cell2mat( ...
        normalizeNumericList(record.values))),1,[]);
elseif state=="COMPLETE"
    error("radia:optuna:ImportStudy", ...
        "Trial %d is COMPLETE but carries no value.",double(record.number));
end
if state=="COMPLETE" && numel(values)~=objectiveCount
    error("radia:optuna:ImportStudy", ...
        "Trial %d has %d values for %d directions.", ...
        double(record.number),numel(values),objectiveCount);
end

intermediate=radia.optuna.Trial.emptyIntermediateTable();
steps=normalizeList(record,"intermediate_values");
for index=1:numel(steps)
    item=steps{index};
    intermediate(end+1,:)={double(item.step),double(item.value), ...
        datetime("now",TimeZone="local")}; %#ok<AGROW>
end

constraintPresent=isfield(record,"constraint_present") && ...
    logical(record.constraint_present);
constraints=zeros(1,0);
if constraintPresent && isfield(record,"constraints") && ...
        ~isempty(record.constraints)
    constraints=reshape(double(cell2mat( ...
        normalizeNumericList(record.constraints))),1,[]);
end
[userAttrs,userAttrNames]=attributeStruct(record,"user_attrs");
[systemAttrs,systemAttrNames]=attributeStruct(record,"system_attrs");
frozen=radia.optuna.FrozenTrial( ...
    Number=double(record.number),State=state,Values=values, ...
    Params=params,Distributions=distributions, ...
    IntermediateValues=intermediate, ...
    UserAttrs=userAttrs,SystemAttrs=systemAttrs, ...
    ConstraintPresent=constraintPresent,Constraints=constraints, ...
    DatetimeStart=parseTimestamp(record,"datetime_start"), ...
    DatetimeComplete=parseTimestamp(record,"datetime_complete"));
end

function [attrs,originalNames]=attributeStruct(record,field)
attrs=struct();
originalNames=struct();
entries=normalizeList(record,field);
if isempty(entries)
    return
end
names=cellfun(@(item)string(item.name),entries);
keys=radia.optuna.Trial.claimKeys(names);
for index=1:numel(entries)
    key=char(keys(index));
    attrs.(key)=jsondecode(char(entries{index}.value_json));
    originalNames.(key)=char(string(entries{index}.name));
end
end

function applyStudyAttributes(study,payload,field,setter)
for entry=reshape(normalizeList(payload,field),1,[])
    item=entry{1};
    setter(string(item.name),jsondecode(char(item.value_json)));
end
end

function items=normalizeList(container,field)
items={};
if ~isfield(container,field) || isempty(container.(field))
    return
end
value=container.(field);
if iscell(value)
    items=reshape(value,1,[]);
elseif isstruct(value)
    items=num2cell(reshape(value,1,[]));
else
    items={value};
end
end

function items=normalizeNumericList(value)
if iscell(value)
    items=reshape(value,1,[]);
else
    items=num2cell(reshape(double(value),1,[]));
end
end

function value=parseTimestamp(record,field)
value=datetime(NaT,TimeZone="local");
if ~isfield(record,field) || isempty(record.(field))
    return
end
value=datetime(string(record.(field)), ...
    InputFormat="uuuu-MM-dd'T'HH:mm:ss.SSSSSS",TimeZone="local");
end
