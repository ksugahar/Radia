# Standalone GUI candidate validation

Use a dedicated wheel-installed Python 3.12 venv without Radia or external Qt.
Do not repoint shared editable installations to run these tests.

1. Run `tools/check_cubit_standalone.py` with isolated Python (`-I`).
   It tests installation/verification in a fake profile, never the real profile.
2. For the explicitly authorized GUI lane, ensure no user Cubit session runs.
   Run `run_standalone_gui_probe.py <new-output-directory>` with that same venv.
   The probe loads menu code from the installed wheel by file path, not by
   importing Python-3.12 native libraries inside Cubit. It registers six actions,
   captures the window and exports an order-2 sphere. The controller requires
   successful probe AND Cubit exit zero; intermediate gui-result.json alone is
   not acceptance. Repeat with a different output directory for a fresh process.
3. Validate each sphere.vol using the structural, CAD and strict body/outer label
   gate in cubit_mesh_export.smoke_test._validate_exported_vol.

Diagnostic environment switches: CME_GUI_TEST_BASELINE=1 skips the candidate;
CME_GUI_TEST_MENU_ONLY=1 skips mesh export. Neither is full GUI acceptance.

The startup-file registration is tested in a scratch profile. The GUI probe
uses -noinitfile with an explicit candidate path: two passing starts do not
prove persistence of a user's imported WorkflowToolbar. Manual package import,
dialog interaction/cancel paths and persisted-toolbar reopen remain separate
acceptance steps before release. Existing Radia legacy GUI assets are retained
pending adapter consolidation; no new version or publication is implied.
