# Standalone GUI candidate validation

Use a dedicated wheel-installed Python 3.12 venv without Radia or external Qt.
Do not repoint shared editable installations to run these tests.

1. Run `tools/check_cubit_standalone.py` with isolated Python (`-I`).
   It tests installation/verification in a fake profile, never the real profile.
2. For the explicitly authorized GUI lane, ensure no user Cubit session runs.
   Run `run_standalone_gui_probe.py <new-output-directory>` with that same venv.
   The probe loads menu code from the installed wheel by file path, not by
   importing Python-3.12 native libraries inside Cubit. It registers six actions,
   captures the window and exports an order-2 sphere. Set
   `CME_GUI_TEST_DIALOGS=1` to exercise six cancel paths, real Netgen accept +
   check-vol result tables, and a missing-checker error dialog. Only the settings
   file and external viewer are redirected; export/check use the real backend.
   The controller requires
   successful probe AND Cubit exit zero; intermediate gui-result.json alone is
   not acceptance. Repeat with a different output directory for a fresh process.
3. Validate each sphere.vol using the structural, CAD and strict body/outer label
   gate in cubit_mesh_export.smoke_test._validate_exported_vol.

Diagnostic environment switches: CME_GUI_TEST_BASELINE=1 skips the candidate;
CME_GUI_TEST_MENU_ONLY=1 skips mesh export. Neither is full GUI acceptance.

4. Run `run_standalone_toolbar_persistence.py <new-output-directory>` with the
   wheel-only interpreter. This uses the real Custom Toolbar Editor to import
   its generated package, exits, then starts a new Cubit process without an
   explicit toolbar load. It verifies all six restored buttons, clicks each,
   cancels each real dialog, and checks that the imported copy supplies the code.
   The toolbar gets its own dock row so all buttons are visible for the test.
   Native clicks are queued, not called synchronously from a nested Python
   callback: otherwise Cubit's Python importer can deadlock behind the driver.

The startup registration gate uses a scratch profile. Official persistence uses
the user's actual Cubit.ini because Qt ignores APPDATA overrides on Windows.
The controller requires no active Cubit session, saves the exact original file
to its output directory, temporarily clears only the toolbar registration, and
restores the original bytes/hash in finally. It never overwrites the user's
imported toolbar directory. A still-running Cubit blocks restoration and must
be resolved before completion; do not discard the saved preferences.

The legacy Radia menu, startup, installer and compatibility bridges are removed.
`cubit-toolbar-smoke-test` and its embedded probe belong to cubit-mesh-export;
their unit tests live in `packages/cubit-mesh-export/tests`. Re-run
`cubit-plugin-install` from the selected external interpreter to replace old
startup references before running the explicitly scoped GUI acceptance.
No new version, publication,
shared editable change or LAB/100 deployment is implied by these local tests.
