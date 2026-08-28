function payload = export_study(study, options)
%EXPORT_STUDY Serialize a study for handoff to an upstream Optuna storage.
%   PAYLOAD = EXPORT_STUDY(STUDY) returns a struct describing every trial in
%   STUDY using Optuna's own distribution JSON, so the Python side can
%   rebuild it with optuna.distributions.json_to_distribution and
%   Study.add_trials.
%
%   EXPORT_STUDY(STUDY, Path=FILE) also writes the payload to FILE.
%
%   This is the interoperability path for behavior radia.optuna does not
%   implement: MATLAB's table/MAT storage is not an optuna.storages backend,
%   so a study cannot simply be opened from Python. Export it, run the
%   Python-only feature (dashboard, visualization, importance, an
%   integration sampler), and import the result back.
%
%   Names are carried as records rather than JSON object keys because
%   parameter and attribute names need not be valid MATLAB field names; the
%   export must not silently rename them.
%
%   See also radia.optuna.import_study, radia.optuna.distribution_to_json.

arguments
    study (1,1) radia.optuna.Study
    options.Path (1,1) string = ""
end

trials = study.TrialTable;
parameters = study.ParamTable;
objectives = study.ObjectiveTable;
intermediate = study.IntermediateTable;
userAttrs = study.UserAttrTable;
systemAttrs = study.SystemAttrTable;

records = cell(1, height(trials));
for index = 1:height(trials)
    number = trials.TrialNumber(index);
    record = struct( ...
        "number", number, ...
        "state", char(trials.State(index)), ...
        "values", {objectiveValues(study, objectives, trials, index, number)}, ...
        "params", {parameterRecords(parameters, number)}, ...
        "user_attrs", {attributeRecords(userAttrs, number)}, ...
        "system_attrs", {attributeRecords(systemAttrs, number)}, ...
        "intermediate_values", {intermediateRecords(intermediate, number)}, ...
        "datetime_start", isoTimestamp(trials.StartTime(index)), ...
        "datetime_complete", isoTimestamp(trials.EndTime(index)));
    records{index} = record;
end

payload = struct( ...
    "schema", "radia.optuna.study-export.v1", ...
    "study_name", char(study.Name), ...
    "directions", {cellstr(reshape(study.Directions, 1, []))}, ...
    "metric_names", {cellstr(reshape(study.MetricNames, 1, []))}, ...
    "user_attrs", {studyAttributes(study.UserAttrs)}, ...
    "system_attrs", {studyAttributes(study.SystemAttrs)}, ...
    "trial_count", height(trials), ...
    "trials", {records});

if strlength(options.Path) > 0
    fid = fopen(options.Path, "w", "n", "UTF-8");
    if fid < 0
        error("radia:optuna:ExportStudy", ...
            "Cannot open '%s' for writing.", options.Path);
    end
    closer = onCleanup(@() fclose(fid));
    fwrite(fid, jsonencode(payload, PrettyPrint=true), "char");
    clear closer
end
end

function values = objectiveValues(study, objectives, trials, index, number)
rows = objectives.TrialNumber == number;
if any(rows)
    ordered = sortrows(objectives(rows, :), "ObjectiveIndex");
    values = num2cell(reshape(ordered.Value, 1, []));
    return
end
% A PRUNED trial owns no objective row; Study.tell stored the value at its
% deepest reported step in the scalar column, exactly as upstream does.
scalar = trials.Value(index);
if isscalar(study.Directions) && ~isnan(scalar)
    values = {scalar};
else
    values = {};
end
end

function records = parameterRecords(parameters, number)
rows = find(parameters.TrialNumber == number);
records = cell(1, numel(rows));
for index = 1:numel(rows)
    row = rows(index);
    if isfinite(parameters.ValueNumeric(row))
        value = parameters.ValueNumeric(row);
    else
        value = jsondecode(parameters.ValueText(row));
    end
    distribution = radia.optuna.internal.DistributionCodec.decode( ...
        parameters.Kind(row), parameters.Distribution(row));
    records{index} = struct( ...
        "name", char(parameters.Name(row)), ...
        "value", value, ...
        "distribution", char( ...
            radia.optuna.distribution_to_json(distribution)));
end
end

function records = attributeRecords(source, number)
rows = find(source.TrialNumber == number);
records = cell(1, numel(rows));
for index = 1:numel(rows)
    row = rows(index);
    records{index} = struct( ...
        "name", char(source.Name(row)), ...
        "value_json", char(source.ValueJSON(row)));
end
end

function records = intermediateRecords(intermediate, number)
rows = find(intermediate.TrialNumber == number);
records = cell(1, numel(rows));
for index = 1:numel(rows)
    row = rows(index);
    records{index} = struct( ...
        "step", intermediate.Step(row), ...
        "value", intermediate.Value(row));
end
end

function records = studyAttributes(attrs)
names = string(fieldnames(attrs));
records = cell(1, numel(names));
for index = 1:numel(names)
    records{index} = struct( ...
        "name", char(names(index)), ...
        "value_json", char(string(jsonencode(attrs.(names(index))))));
end
end

function text = isoTimestamp(value)
if isnat(value)
    text = [];
    return
end
text = char(string(value, "uuuu-MM-dd'T'HH:mm:ss.SSSSSS"));
end
