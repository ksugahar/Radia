function study = import_study(payload, options)
%IMPORT_STUDY Rebuild a study from a radia.optuna.study-export.v1 payload.
%   STUDY = IMPORT_STUDY(PATH) reads the JSON document written by
%   radia.optuna.export_study, or by the Python bridge after upstream Optuna
%   has worked on it, and returns a study holding the same trials.
%
%   STUDY = IMPORT_STUDY(PAYLOAD) accepts an already-decoded struct.
%
%   Trial numbers are reassigned in ascending order of the exported numbers,
%   which is what Study.add_trial does; the payload's own numbers are kept in
%   each trial's system attributes so a round trip stays traceable.
%
%   See also radia.optuna.export_study.

arguments
    payload
    options.Name (1,1) string = ""
    options.StoragePath (1,1) string = ""
    options.Sampler = []
    options.AutoSave (1,1) logical = true
end

if isstring(payload) || ischar(payload)
    path = string(payload);
    if ~isfile(path)
        error("radia:optuna:ImportStudy", ...
            "Export file '%s' does not exist.", path);
    end
    payload = jsondecode(fileread(path));
end
if ~isstruct(payload) || ~isscalar(payload)
    error("radia:optuna:ImportStudy", ...
        "An export payload must be a scalar struct or a file path.");
end
if ~isfield(payload, "schema") || ...
        string(payload.schema) ~= "radia.optuna.study-export.v1"
    error("radia:optuna:ImportStudy", ...
        "Unsupported export schema; expected radia.optuna.study-export.v1.");
end

name = options.Name;
if strlength(name) == 0 && isfield(payload, "study_name")
    name = string(payload.study_name);
end
directions = reshape(string(payload.directions), 1, []);
study = radia.optuna.Study(Name=name, Directions=directions, ...
    Sampler=options.Sampler, StoragePath=options.StoragePath, ...
    AutoSave=options.AutoSave);

if isfield(payload, "metric_names") && ~isempty(payload.metric_names)
    study.set_metric_names(reshape(string(payload.metric_names), 1, []));
end
applyStudyAttributes(study, payload, "user_attrs", @study.set_user_attr);
applyStudyAttributes(study, payload, "system_attrs", @study.set_system_attr);

trials = normalizeList(payload, "trials");
numbers = cellfun(@(t) double(t.number), trials);
[~, order] = sort(numbers);
for index = reshape(order, 1, [])
    [frozen, originalNames] = frozenFromRecord(trials{index}, numel(directions));
    study.addTrial(frozen, OriginalNames=originalNames);
end
end

function [frozen, originalNames] = frozenFromRecord(record, objectiveCount)
params = struct();
distributions = struct();
originalNames = struct();
% A FrozenTrial can only hold MATLAB-valid field names, so remember the
% name each key stands for and hand it to addTrial; otherwise "x-1" would
% come back as "x_1" and the export would not round-trip.
claimed = radia.optuna.Trial.claimKeys( ...
    cellfun(@(item) string(item.name), normalizeList(record, "params")));
entries = normalizeList(record, "params");
for index = 1:numel(entries)
    item = entries{index};
    key = char(claimed(index));
    params.(key) = item.value;
    distributions.(key) = radia.optuna.json_to_distribution( ...
        string(item.distribution));
    originalNames.(key) = char(string(item.name));
end

state = string(record.state);
values = NaN;
if isfield(record, "values") && ~isempty(record.values)
    values = reshape(double(cell2mat(normalizeNumericList(record.values))), 1, []);
elseif state == "COMPLETE"
    error("radia:optuna:ImportStudy", ...
        "Trial %d is COMPLETE but carries no value.", double(record.number));
end
if state == "COMPLETE" && numel(values) ~= objectiveCount
    error("radia:optuna:ImportStudy", ...
        "Trial %d has %d values for %d directions.", ...
        double(record.number), numel(values), objectiveCount);
end

intermediate = radia.optuna.Trial.emptyIntermediateTable();
steps = normalizeList(record, "intermediate_values");
for index = 1:numel(steps)
    item = steps{index};
    intermediate(end+1,:) = {double(item.step), double(item.value), ...
        datetime("now", TimeZone="local")}; %#ok<AGROW>
end

arguments_ = { ...
    "State", state, ...
    "Values", values, ...
    "Params", params, ...
    "Distributions", distributions, ...
    "IntermediateValues", intermediate, ...
    "UserAttrs", attributeStruct(record, "user_attrs"), ...
    "SystemAttrs", attributeStruct(record, "system_attrs")};
frozen = radia.optuna.FrozenTrial(arguments_{:}, ...
    Number=double(record.number), ...
    DatetimeStart=parseTimestamp(record, "datetime_start"), ...
    DatetimeComplete=parseTimestamp(record, "datetime_complete"));
end

function attrs = attributeStruct(record, field)
attrs = struct();
for entry = reshape(normalizeList(record, field), 1, [])
    item = entry{1};
    attrs.(matlab.lang.makeValidName(string(item.name))) = ...
        jsondecode(char(item.value_json));
end
end

function applyStudyAttributes(study, payload, field, setter)
for entry = reshape(normalizeList(payload, field), 1, [])
    item = entry{1};
    setter(string(item.name), jsondecode(char(item.value_json)));
end
end

function items = normalizeList(container, field)
% jsondecode collapses a homogeneous array of objects into a struct array
% and an empty array into [], so restore a cell list either way.
items = {};
if ~isfield(container, field) || isempty(container.(field))
    return
end
value = container.(field);
if iscell(value)
    items = reshape(value, 1, []);
elseif isstruct(value)
    items = num2cell(reshape(value, 1, []));
else
    items = {value};
end
end

function items = normalizeNumericList(value)
if iscell(value)
    items = reshape(value, 1, []);
else
    items = num2cell(reshape(double(value), 1, []));
end
end

function value = parseTimestamp(record, field)
value = datetime(NaT, TimeZone="local");
if ~isfield(record, field) || isempty(record.(field))
    return
end
value = datetime(string(record.(field)), ...
    InputFormat="uuuu-MM-dd'T'HH:mm:ss.SSSSSS", TimeZone="local");
end
