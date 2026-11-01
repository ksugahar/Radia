function tests = test_optuna_hdf5_storage
%TEST_OPTUNA_HDF5_STORAGE Typed HDF5 Study/session persistence contracts.
%   Guards: Study and OptimizationSession persist only typed HDF5 (no MAT
%   path, reader or fallback); every StudyData field roundtrips with exact
%   types, order, UTF-8 text, NaN/Inf and datetimes; the sampler RNG state
%   (uint32[624], index, has_gauss, gauss bits) survives save/load so a
%   resumed TPE study proposes exactly what an uninterrupted one does;
%   schema mismatches and non-HDF5 inputs fail loudly; session callbacks
%   are recorded only as unresolved descriptions and are never restored.
%   Assumes: MATLAB HDF5 low-level API; no tolerance is involved.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% An installed-wheel run names its MATLAB directory explicitly; otherwise
% the source tree beside this test is exercised.
matlabDirectory = string(getenv("RADIA_OPTUNA_MATLAB_PATH"));
if matlabDirectory == ""
    root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    matlabDirectory = fullfile(root, "packages", "radia-optuna", "matlab");
end
entries = string(strsplit(path, pathsep));
testCase.TestData.RemovePath = ~any(strcmpi(entries, string(matlabDirectory)));
if testCase.TestData.RemovePath
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testStudyTypedRoundTrip(testCase)
path = temporaryPath(testCase);
text = string(char([26085 26412 35486 32 955]));
study = radia.optuna.Study(Name="typed " + text, StoragePath=path, ...
    Sampler=radia.optuna.RandomSampler(5));
study.setUserAttr("label", text);
study.setMetricNames("loss");
study.enqueueTrial(struct("x", 0.25, "k", 2, "c", "b"));
for index = 1:4
    trial = study.ask();
    trial.suggestFloat("x", -1, 1);
    trial.suggestInteger("k", 0, 4);
    trial.suggestCategorical("c", ["a", "b", "c"]);
    trial.setUserAttr("note", text + index);
    trial.report(index, 0);
    switch index
        case 2
            study.fail(trial, "failed " + text);
        case 3
            study.tell(trial, Inf);
        otherwise
            study.tell(trial, index - 0.5);
    end
end
running = study.ask();
running.suggestFloat("x", -1, 1);
study.save();

restored = radia.optuna.Study(StoragePath=path, AutoSave=false);
verifyEqual(testCase, restored.Name, study.Name);
verifyEqual(testCase, restored.MetricNames, study.MetricNames);
verifyEqual(testCase, restored.UserAttrs, study.UserAttrs);
tables = ["TrialTable","ParamTable","IntermediateTable","UserAttrTable", ...
    "SystemAttrTable","ConstraintTable","ConstraintCountTable", ...
    "ObjectiveTable","SamplerStateTable","QueueParamTable"];
for name = tables
    verifyTrue(testCase, isequaln(restored.(name), study.(name)), ...
        "Table " + name + " changed in the HDF5 roundtrip.");
end
verifyEqual(testCase, restored.TrialTable.ErrorMessage(2), "failed " + text);
verifyEqual(testCase, restored.TrialTable.Value(3), Inf);
verifyEqual(testCase, restored.TrialTable.State(end), "RUNNING");

% Format: real HDF5 with explicit schema, units and index conventions,
% read here through the high-level API independent of the writer.
verifyTrue(testCase, logical(H5F.is_hdf5(char(path))));
verifyEqual(testCase, string(h5readatt(path, "/", "schema")), ...
    "radia.optuna.study");
verifyEqual(testCase, h5readatt(path, "/", "schema_version"), int64(7));
verifyEqual(testCase, string(h5readatt(path, "/", "encoding")), ...
    "radia.optuna.typed-tree.v1");
verifyEqual(testCase, string(h5readatt(path, "/StudyData", ...
    "units:TrialTable.Duration_s")), "s");
verifyEqual(testCase, h5readatt(path, "/StudyData", ...
    "index_base:TrialTable.TrialNumber"), int64(0));
verifyEqual(testCase, h5readatt(path, "/StudyData", ...
    "index_base:ObjectiveTable.ObjectiveIndex"), int64(1));
verifyClass(testCase, h5read(path, "/StudyData/uint32"), "uint32");
verifyClass(testCase, h5read(path, "/StudyData/text_bytes"), "uint8");
fields = radia.optuna.Study.readStorageFile(path);
verifyEqual(testCase, fields.NextTrialNumber, 5);
verifyEqual(testCase, string(fieldnames(fields)).', ["Schema","Version", ...
    "Name","Directions","NextTrialNumber","TrialTable","ParamTable", ...
    "IntermediateTable","TrialIntermediateSnapshots","UserAttrTable", ...
    "SystemAttrTable","ConstraintTable","ConstraintCountTable", ...
    "UserAttrs","SystemAttrs","MetricNames","ObjectiveTable", ...
    "SamplerStateTable","QueueParamTable"]);
end

function testResumedTpeMatchesUninterruptedRun(testCase)
path = temporaryPath(testCase);
objective = @(trial) objectiveValue(trial);
% The uninterrupted reference is persisted too, so it records the same
% sampler-state rows that the resumed study restores from HDF5.
reference = radia.optuna.Study(StoragePath=temporaryPath(testCase), ...
    Sampler=radia.optuna.TPESampler(Seed=314, NStartupTrials=3));
reference.optimize(objective, 10);

first = radia.optuna.Study(StoragePath=path, Sampler= ...
    radia.optuna.TPESampler(Seed=314, NStartupTrials=3));
first.optimize(objective, 5);
stored = radia.optuna.Study.readStorageFile(path);
memoryStates = first.SamplerStateTable.State;
fileStates = stored.SamplerStateTable.State;
memoryState = memoryStates{end}.random_state;
fileState = fileStates{end}.random_state;
verifyClass(testCase, fileState.mt, "uint32");
verifyEqual(testCase, numel(fileState.mt), 624);
verifyEqual(testCase, fileState.mt, memoryState.mt);
verifyEqual(testCase, fileState.index, memoryState.index);
verifyEqual(testCase, fileState.has_gauss, memoryState.has_gauss);
verifyEqual(testCase, typecast(fileState.gauss, "uint64"), ...
    typecast(memoryState.gauss, "uint64"));
clear first

resumed = radia.optuna.Study(StoragePath=path, Sampler= ...
    radia.optuna.TPESampler(Seed=314, NStartupTrials=3));
resumed.optimize(objective, 5);
verifyEqual(testCase, resumed.ParamTable.ValueNumeric, ...
    reference.ParamTable.ValueNumeric);
verifyEqual(testCase, resumed.TrialTable.Value, reference.TrialTable.Value);
verifyTrue(testCase, isequaln(resumed.SamplerStateTable.State, ...
    reference.SamplerStateTable.State));
end

function testSchemaAndFormatMismatchFailLoudly(testCase)
path = temporaryPath(testCase);
study = radia.optuna.Study(StoragePath=path, ...
    Sampler=radia.optuna.RandomSampler(1));
trial = study.ask();
trial.suggestFloat("x", 0, 1);
study.tell(trial, 1);
verifyError(testCase, @() radia.optuna.Study( ...
    StoragePath=string(tempname("C:\temp")) + ".mat"), ...
    "radia:optuna:StorageFormat");

h5writeatt(path, "/", "schema_version", int64(6));
verifyError(testCase, @() radia.optuna.Study.readStorageFile(path), ...
    "radia:optuna:StorageVersion");
copyfile(path + ".bak", path, "f");
h5writeatt(path, "/", "schema", "radia.optuna.other");
verifyError(testCase, @() radia.optuna.Study.readStorageFile(path), ...
    "radia:optuna:StorageVersion");

notHdf5 = temporaryPath(testCase);
StudyData = struct("Schema", "radia.optuna.study"); %#ok<NASGU>
save(notHdf5, "StudyData", "-mat");
verifyError(testCase, @() radia.optuna.Study.readStorageFile(notHdf5), ...
    "radia:optuna:StorageFormat");
verifyError(testCase, @() radia.optuna.Study(StoragePath=notHdf5, ...
    AutoSave=false), "radia:optuna:StorageRecovery");

verifyError(testCase, @() radia.optuna.internal.encodeTypedTree( ...
    struct("callback", {{@sin}})), "radia:optuna:StorageUnsupported");
end

function testSessionCallbacksStayUnresolved(testCase)
storage = temporaryPath(testCase);
parameters = radia.optuna.OptimizationParameter("x", Minimum=0, Maximum=1);
parameters(2) = radia.optuna.OptimizationParameter("c", ...
    Type="categorical", Choices={"left", 2});
options = radia.optuna.optimoptions(MaxTrials=3, Sampler="random", ...
    Seed=11, Display="off", StoragePath=storage, MaxTime=Inf);
session = radia.optuna.OptimizationSession(@(v) v.x, parameters, options);
session.start();
session.runNext();
session.pause();
sessionPath = session.save();
addTeardown(testCase, @() deleteIfPresent(sessionPath));
verifyEqual(testCase, sessionPath, ...
    radia.optuna.OptimizationSession.sessionPathFor(storage));
verifyTrue(testCase, endsWith(sessionPath, ".session.h5"));

data = radia.optuna.OptimizationSession.readSessionFile(sessionPath);
verifyEqual(testCase, data.Parameters, session.Parameters);
verifyEqual(testCase, data.Options, session.Options);
verifyEqual(testCase, data.UnresolvedCallbacks.Slot, "ObjectiveFcn");
verifyEqual(testCase, data.UnresolvedCallbacks.Description, "@(v)v.x");
verifyEqual(testCase, h5readatt(sessionPath, "/SessionData", ...
    "restorable:UnresolvedCallbacks"), int64(0));

verifyError(testCase, @() radia.optuna.OptimizationSession.load( ...
    sessionPath), "radia:optuna:SessionCallbackUnresolved");
restored = radia.optuna.OptimizationSession.load(sessionPath, ...
    ObjectiveFcn=@(v) v.x);
verifyEqual(testCase, restored.State, "paused");
verifyEqual(testCase, restored.snapshot().complete, 1);

callbackStorage = temporaryPath(testCase);
callbackOptions = radia.optuna.optimoptions(MaxTrials=2, ...
    Sampler="random", Seed=3, Display="off", ...
    StoragePath=callbackStorage, OutputFcn=@(varargin) false);
withCallback = radia.optuna.OptimizationSession(@(v) v.x, ...
    parameters(1), callbackOptions);
callbackPath = withCallback.save();
addTeardown(testCase, @() deleteIfPresent(callbackPath));
recorded = radia.optuna.OptimizationSession.readSessionFile(callbackPath);
verifyEqual(testCase, recorded.UnresolvedCallbacks.Slot, ...
    ["ObjectiveFcn"; "OutputFcn"]);
verifyError(testCase, @() radia.optuna.OptimizationSession.load( ...
    callbackPath, ObjectiveFcn=@(v) v.x), ...
    "radia:optuna:SessionCallbackUnresolved");
end

function testSessionRebindsAllCallbacksAndResumesIdentically(testCase)
parameter = radia.optuna.OptimizationParameter("x", Minimum=0, Maximum=1);
objective = @(v) (v.x - 0.3)^2;
output = @(varargin) false;
plotA = @(varargin) false;
plotB = @(varargin) false;
constraint = @(v) v.x - 0.8;
optionsFor = @(storage) radia.optuna.optimoptions(MaxTrials=5, ...
    Sampler="random", Seed=41, Display="off", StoragePath=storage, ...
    OutputFcn=output, PlotFcn={plotA, plotB}, ConstraintFcn=constraint);

reference = radia.optuna.OptimizationSession(objective, parameter, ...
    optionsFor(temporaryPath(testCase)));
addTeardown(testCase, @() deleteIfPresent( ...
    radia.optuna.OptimizationSession.sessionPathFor(reference.Options.StoragePath)));
reference.start();
reference.run();

storage = temporaryPath(testCase);
first = radia.optuna.OptimizationSession(objective, parameter, ...
    optionsFor(storage));
first.start();
first.runNext();
first.runNext();
first.pause();
sessionPath = first.save();
addTeardown(testCase, @() deleteIfPresent(sessionPath));
clear first
recorded = radia.optuna.OptimizationSession.readSessionFile(sessionPath);
verifyEqual(testCase, recorded.UnresolvedCallbacks.Slot, ...
    ["ObjectiveFcn"; "OutputFcn"; "PlotFcn"; "PlotFcn"; "ConstraintFcn"]);
verifyEqual(testCase, recorded.UnresolvedCallbacks.Position, [1; 1; 1; 2; 1]);

resumed = radia.optuna.OptimizationSession.load(sessionPath, ...
    ObjectiveFcn=objective, OutputFcn=output, PlotFcn={plotA, plotB}, ...
    ConstraintFcn=constraint);
verifyEqual(testCase, resumed.Options.OutputFcn, {output});
verifyEqual(testCase, resumed.Options.PlotFcn, {plotA, plotB});
verifyEqual(testCase, resumed.Options.ConstraintFcn, constraint);
resumed.resume();
resumed.run();
verifyEqual(testCase, resumed.State, "completed");
verifyEqual(testCase, resumed.Study.ParamTable.ValueNumeric, ...
    reference.Study.ParamTable.ValueNumeric);
verifyEqual(testCase, resumed.Study.TrialTable.Value, ...
    reference.Study.TrialTable.Value);
verifyEqual(testCase, resumed.Study.ConstraintTable.Value, ...
    reference.Study.ConstraintTable.Value);
end

function testSessionRebindRejectsMissingWrongOrExtraCallbacks(testCase)
parameter = radia.optuna.OptimizationParameter("x", Minimum=0, Maximum=1);
objective = @(v) v.x;
output = @(varargin) false;
storage = temporaryPath(testCase);
options = radia.optuna.optimoptions(MaxTrials=2, Sampler="random", ...
    Seed=5, Display="off", StoragePath=storage, OutputFcn=output);
session = radia.optuna.OptimizationSession(objective, parameter, options);
sessionPath = session.save();
addTeardown(testCase, @() deleteIfPresent(sessionPath));
id = "radia:optuna:SessionCallbackUnresolved";
load = @(varargin) radia.optuna.OptimizationSession.load( ...
    sessionPath, varargin{:});
% Missing required callbacks.
verifyError(testCase, @() load(OutputFcn=output), id);
verifyError(testCase, @() load(ObjectiveFcn=objective), id);
% Wrong count for a recorded slot.
verifyError(testCase, @() load(ObjectiveFcn=objective, ...
    OutputFcn={output, output}), id);
% Callback supplied for a slot that was not recorded.
verifyError(testCase, @() load(ObjectiveFcn=objective, ...
    OutputFcn=output, PlotFcn=output), id);
verifyError(testCase, @() load(ObjectiveFcn=objective, ...
    OutputFcn=output, ConstraintFcn=@(v) v.x), id);
% Wrong type is not accepted as a callback.
verifyError(testCase, @() load(ObjectiveFcn=objective, ...
    OutputFcn="@(varargin)false"), id);
restored = load(ObjectiveFcn=objective, OutputFcn=output);
verifyEqual(testCase, restored.Options.OutputFcn, {output});
end

function testConstructorResumeEnforcesRecordedCallbacks(testCase)
parameter = radia.optuna.OptimizationParameter("x", Minimum=0, Maximum=1);
objective = @(v) v.x;
output = @(varargin) false;
storage = temporaryPath(testCase);
sessionPath = radia.optuna.OptimizationSession.sessionPathFor(storage);
addTeardown(testCase, @() deleteIfPresent(sessionPath));
first = radia.optuna.OptimizationSession(objective, parameter, ...
    radia.optuna.optimoptions(MaxTrials=2, Sampler="random", Seed=7, ...
    Display="off", StoragePath=storage, OutputFcn=output));
first.start();
first.runNext();
first.pause();
clear first
before = fileread_bytes(sessionPath);

withoutOutput = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=7, Display="off", StoragePath=storage);
verifyError(testCase, @() radia.optuna.OptimizationSession(objective, ...
    parameter, withoutOutput), "radia:optuna:SessionCallbackUnresolved");
verifyEqual(testCase, fileread_bytes(sessionPath), before, ...
    "A rejected constructor resume must not rewrite the checkpoint.");
extraPlot = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=7, Display="off", StoragePath=storage, OutputFcn=output, ...
    PlotFcn=output);
verifyError(testCase, @() radia.optuna.OptimizationSession(objective, ...
    parameter, extraPlot), "radia:optuna:SessionCallbackUnresolved");

extended = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=7, Display="off", StoragePath=storage, OutputFcn=output);
resumed = radia.optuna.OptimizationSession(objective, parameter, extended);
verifyEqual(testCase, resumed.Options.MaxTrials, 4);
verifyEqual(testCase, resumed.Options.OutputFcn, {output});
verifyEqual(testCase, resumed.snapshot().complete, 1);
recorded = radia.optuna.OptimizationSession.readSessionFile(sessionPath);
verifyEqual(testCase, recorded.UnresolvedCallbacks.Slot, ...
    ["ObjectiveFcn"; "OutputFcn"]);
verifyEqual(testCase, recorded.Options.MaxTrials, 4);
end

function testUnsupportedDiagnosticNamesLabelNotSearchPath(testCase)
try
    radia.optuna.internal.encodeTypedTree(struct("callback", {{@sin}}), ...
        "ProbeLabel");
    verifyFail(testCase, "encodeTypedTree accepted a function handle.");
catch failure
    verifyEqual(testCase, string(failure.identifier), ...
        "radia:optuna:StorageUnsupported");
    message = string(failure.message);
    verifyTrue(testCase, contains(message, "ProbeLabel"), message);
    verifyTrue(testCase, contains(message, "node "), message);
    verifyFalse(testCase, contains(message, pathsep) || ...
        contains(message, matlabroot) || contains(message, ":\"), message);
end
end

function bytes = fileread_bytes(path)
handle = fopen(path, "r");
closer = onCleanup(@() fclose(handle));
bytes = fread(handle, Inf, "*uint8");
end

function value = objectiveValue(trial)
x = trial.suggestFloat("x", -2, 2);
y = trial.suggestFloat("y", -2, 2);
value = (x - 0.3)^2 + (y + 0.7)^2;
end

function path = temporaryPath(testCase)
path = string(tempname("C:\temp")) + ".h5";
addTeardown(testCase, @() deleteIfPresent(path));
addTeardown(testCase, @() deleteIfPresent(path + ".bak"));
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
