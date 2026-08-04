function tests = test_ih_geometry_update
%TEST_IH_GEOMETRY_UPDATE Geometry watch/rebuild block and engine.
%   Mechanics under test: content-fingerprint staleness detection, the
%   assemble command running exactly when needed, configuration reload
%   into the model workspace, crossed .vol/.step repair, and the
%   fail-loud paths (stale with auto off, two STEP inputs).
tests = functiontests(localfunctions);
end

function [model, closer] = freshModel()
model = "geom_upd_" + erase(string(java.util.UUID.randomUUID), "-");
new_system(model);
closer = onCleanup(@() close_system(model, 0));
end

function [work, cleanupDir, wp, coil, golden, configFile, runsFile, command] = fixture()
work = string(tempname("C:\temp"));
mkdir(work);
cleanupDir = onCleanup(@() rmdir(work, "s"));
wp = fullfile(work, "wp.vol");
writelines("mesh v1", wp);
coil = fullfile(work, "coil.step");
writelines("step v1", coil);
golden = fullfile(work, "golden.mat");
config = radia.simulink.makeIHNativeSmokeConfig();
save(golden, "config");
configFile = fullfile(work, "config.mat");
runsFile = fullfile(work, "runs.txt");
command = sprintf('copy /y "%s" "%s" >nul && echo run>> "%s"', ...
    golden, configFile, runsFile);
end

function count = runCount(runsFile)
count = 0;
if isfile(runsFile)
    count = numel(readlines(runsFile)) - 1;  % trailing newline
end
end

function configureBlock(model, wp, coil, command, configFile, auto)
block = model + "/Geometry Update";
set_param(block, "wp_vol", char(wp), "coil_file", char(coil), ...
    "assemble_command", char(command), "config_file", char(configFile), ...
    "auto_rebuild", auto);
end

function testRebuildOnFirstRunSkipWhenFreshRebuildOnChange(testCase)
[~, cleanupDir, wp, coil, ~, configFile, runsFile, command] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
configureBlock(model, wp, coil, command, configFile, "on");

first = radia.simulink.updateIHGeometry(model);
verifyTrue(testCase, first.rebuilt);
verifyEqual(testCase, first.revision, 1);
verifyTrue(testCase, isfile(configFile));
verifyTrue(testCase, isfile(configFile + ".fingerprint.json"));
verifyEqual(testCase, runCount(runsFile), 1);
workspace = get_param(model, "ModelWorkspace");
verifyTrue(testCase, workspace.hasVariable("radia_ih_config"));

second = radia.simulink.updateIHGeometry(model);
verifyFalse(testCase, second.rebuilt);
verifyEqual(testCase, string(second.reason), "up-to-date");
verifyEqual(testCase, runCount(runsFile), 1);

writelines("mesh v2 -- user edited", wp);
third = radia.simulink.updateIHGeometry(model);
verifyTrue(testCase, third.rebuilt);
verifyEqual(testCase, third.revision, 2);
verifyEqual(testCase, runCount(runsFile), 2);
end

function testCrossedBoxesRepairedAndInitFcnInstalled(testCase)
[~, cleanupDir, wp, coil, ~, configFile, ~, command] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
verifyTrue(testCase, contains( ...
    string(get_param(model, "InitFcn")), ...
    "radia.simulink.updateIHGeometry"));
% The coil STEP typed into the workpiece box and vice versa.
configureBlock(model, coil, wp, command, configFile, "on");
warned = warning("off", "radia:simulink:IHGeometryRolesReassigned");
restore = onCleanup(@() warning(warned));
status = radia.simulink.updateIHGeometry(model);
verifyTrue(testCase, status.rebuilt);
verifyEqual(testCase, string(status.files(1)), string(wp));
verifyEqual(testCase, string(status.files(2)), string(coil));
verifyNotEmpty(testCase, status.notes);
end

function testBrowseButtonsAreRegisteredBesideGeometryFields(testCase)
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);
mask = Simulink.Mask.get(block);
wpButton = mask.getDialogControl("browse_wp_vol");
coilButton = mask.getDialogControl("browse_coil_file");
verifyNotEmpty(testCase, wpButton);
verifyNotEmpty(testCase, coilButton);
verifyEqual(testCase, string(wpButton.Prompt), "Browse...");
verifyEqual(testCase, string(coilButton.Prompt), "Browse...");
verifyEqual(testCase, string(wpButton.Row), "current");
verifyEqual(testCase, string(coilButton.Row), "current");
verifyTrue(testCase, contains(string(wpButton.Callback), "workpiece"));
verifyTrue(testCase, contains(string(coilButton.Callback), "coil"));
end

function testBrowseAssignsAbsoluteWorkpieceAndCoilPaths(testCase)
[~, cleanupDir, wp, coil] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);

selectedWp = radia.simulink.browseIHGeometryFile( ...
    block, "workpiece", DialogFcn=@(~, ~, ~) selectPath(wp));
selectedCoil = radia.simulink.browseIHGeometryFile( ...
    block, "coil", DialogFcn=@(~, ~, ~) selectPath(coil));

verifyEqual(testCase, string(get_param(block, "wp_vol")), selectedWp);
verifyEqual(testCase, string(get_param(block, "coil_file")), selectedCoil);
verifyTrue(testCase, isfile(selectedWp));
verifyTrue(testCase, isfile(selectedCoil));
verifyTrue(testCase, isAbsolutePath(selectedWp));
verifyTrue(testCase, isAbsolutePath(selectedCoil));
end

function testBrowseCancelLeavesPathUnchanged(testCase)
[~, cleanupDir, wp] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);
set_param(block, "wp_vol", char(wp));
selected = radia.simulink.browseIHGeometryFile( ...
    block, "workpiece", DialogFcn=@(~, ~, ~) deal(0, 0));
verifyEqual(testCase, selected, "");
verifyEqual(testCase, string(get_param(block, "wp_vol")), string(wp));
end

function testBrowseAcceptsCompressedWorkpieceAndMeshedCoil(testCase)
[work, cleanupDir] = fixture(); %#ok<ASGLU>
wp = fullfile(work, "wp.vol.gz");
coil = fullfile(work, "coil.vol");
writelines("compressed mesh placeholder", wp);
writelines("coil mesh", coil);
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);

selectedWp = radia.simulink.browseIHGeometryFile( ...
    block, "workpiece", DialogFcn=@(~, ~, ~) selectPath(wp));
selectedCoil = radia.simulink.browseIHGeometryFile( ...
    block, "coil", DialogFcn=@(~, ~, ~) selectPath(coil));

verifyTrue(testCase, endsWith(lower(selectedWp), ".vol.gz"));
verifyTrue(testCase, endsWith(lower(selectedCoil), ".vol"));
verifyEqual(testCase, string(get_param(block, "wp_vol")), selectedWp);
verifyEqual(testCase, string(get_param(block, "coil_file")), selectedCoil);
end

function testBrowseRejectsWrongExtension(testCase)
[work, cleanupDir] = fixture(); %#ok<ASGLU>
wrong = fullfile(work, "workpiece.msh");
writelines("mesh", wrong);
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);
verifyError(testCase, @() radia.simulink.browseIHGeometryFile( ...
    block, "workpiece", DialogFcn=@(~, ~, ~) selectPath(wrong)), ...
    "radia:simulink:IHGeometryBrowseExtension");
end

function testBrowseRejectsMissingFile(testCase)
[work, cleanupDir] = fixture(); %#ok<ASGLU>
missing = fullfile(work, "missing.step");
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);
verifyError(testCase, @() radia.simulink.browseIHGeometryFile( ...
    block, "coil", DialogFcn=@(~, ~, ~) selectPath(missing)), ...
    "radia:simulink:IHGeometryBrowseMissing");
end

function testBrowseRejectsUnknownRole(testCase)
[model, closer] = freshModel(); %#ok<ASGLU>
block = radia.simulink.addIHGeometryUpdateBlock(model);
verifyError(testCase, ...
    @() radia.simulink.browseIHGeometryFile(block, "rotor"), ...
    "radia:simulink:IHGeometryBrowseRole");
end

function testStaleWithAutoOffErrors(testCase)
[~, cleanupDir, wp, coil, ~, configFile, ~, command] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
configureBlock(model, wp, coil, command, configFile, "off");
verifyError(testCase, ...
    @() radia.simulink.updateIHGeometry(model), ...
    "radia:simulink:IHGeometryUpdateStale");
end

function testForceRebuildOverridesAutoOff(testCase)
[~, cleanupDir, wp, coil, ~, configFile, runsFile, command] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
configureBlock(model, wp, coil, command, configFile, "off");
status=radia.simulink.updateIHGeometry(model,Force=true);
verifyTrue(testCase,status.rebuilt);
verifyEqual(testCase,status.revision,1);
verifyEqual(testCase,runCount(runsFile),1);
end

function testSameContentAtNewPathRebuildsForProvenance(testCase)
[work, cleanupDir, wp, coil, ~, configFile, runsFile, command] = fixture(); %#ok<ASGLU>
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
configureBlock(model,wp,coil,command,configFile,"on");
radia.simulink.updateIHGeometry(model);
replacement=fullfile(work,"wp_repointed.vol");copyfile(wp,replacement);
configureBlock(model,replacement,coil,command,configFile,"on");
status=radia.simulink.updateIHGeometry(model);
verifyTrue(testCase,status.rebuilt);
verifyEqual(testCase,status.revision,2);
verifyEqual(testCase,runCount(runsFile),2);
end

function testUnconfiguredBlockIsInert(testCase)
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
status = radia.simulink.updateIHGeometry(model);
verifyFalse(testCase, status.engaged);
verifyEqual(testCase, string(status.reason), "unconfigured");
end

function testTwoStepInputsError(testCase)
[work, cleanupDir, ~, coil, ~, configFile, ~, command] = fixture(); %#ok<ASGLU>
other = fullfile(work, "second.step");
writelines("step v1", other);
[model, closer] = freshModel(); %#ok<ASGLU>
radia.simulink.addIHGeometryUpdateBlock(model);
configureBlock(model, other, coil, command, configFile, "on");
verifyError(testCase, ...
    @() radia.simulink.updateIHGeometry(model), ...
    "radia:simulink:IHGeometryUpdateTwoSteps");
end

function [fileName, folderName] = selectPath(path)
[folderName, stem, extension] = fileparts(path);
fileName = stem + extension;
folderName = folderName + filesep;
end

function tf = isAbsolutePath(path)
tf = ~isempty(regexp(char(path), '^[A-Za-z]:[\\/]', 'once')) || ...
    startsWith(path, "\\");
end
