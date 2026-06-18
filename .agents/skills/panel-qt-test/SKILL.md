---
name: panel-qt-test
description: Run / write headless PySide6 GUI tests for Radia panels (IH, accel, etc.). Use AFTER editing any radia_*.py panel or its calc_*.py subprocess script to catch behaviour regressions that string-grep tests miss (empty Method combo, hidden widgets feeding build_command, Open GMSH grey out, save/restore by index, calc_*.py argparse rejecting GUI commands, etc.).
---

# panel-qt-test

Headless PySide6 GUI tests for Radia panels. The infrastructure
lives in `tests/panels/` and runs on the Qt `offscreen` platform
plugin so it works on a CI runner without a display.

## When to use this skill

After ANY edit to:
- `src/radia/radia_*.py`           (IHPanel, AccelPanel, EMPanel, ...)
- `src/radia/radia_gui_base.py`    (ModePanel, AnalysisWindow shared)
- `src/radia/panels/calc_*.py`     (subprocess scripts)
- `src/radia/panels/calc_common.py`/`calc_heating_bem.py` (shared helpers)

Especially BEFORE deploying to 100号機 / mdx — the deploy skill
runs the same checks but locally is 5 seconds vs the 30-second
SMB copy + remote python verify roundtrip.

## What it catches

Pitfalls from `panel_gui_pitfalls` MCP knowledge that the existing
string-grep tests in `tests/panels/test_panel_ui_logic.py` do NOT
detect:

| Pitfall | Test that catches it |
|---|---|
| empty Method combo on first launch (`combo_state`) | `test_default_method_is_BEM`, `test_legacy_index_out_of_range_keeps_default` |
| hidden widget feeding `build_command` (`mode_switch`) | `test_FEM_workpiece_widgets_visible` |
| GUI sends `--material custom` to a parser that does not accept it (`subprocess_args`) | `test_FEM[SIBC-sibc]` (command roundtrip via real argparse) |
| Open GMSH grey out after FEM (`silent_except`, `result_keys`) | `test_msh_file_enables`, `test_empty_string_does_not_enable` |
| Mode switch leaks BEM workpiece state into FEM | `test_FEM_solver_items`, `test_BEM_shows_BEM_only_widgets` |
| Run button stays disabled after Browse... (`vol_path` re-inspect) | `test_browse_to_valid_vol_enables_run`, `test_init_uses_restored_vol_path_not_constructor_arg`, `test_browse_vol_invokes_hook` (kubota 2026-05-08, v4.28.1) |

## Run the tests

```bash
cd s:/Radia/01_GitHub
QT_QPA_PLATFORM=offscreen python -m pytest tests/panels/test_ih_panel_qt.py \
                                            tests/panels/test_panel_state_restore.py \
                                            tests/panels/test_open_gmsh_button.py \
                                            tests/panels/test_run_button_browse.py \
                                            -v
```

Expected: **40+ passed in ~5 s**.

The `QT_QPA_PLATFORM=offscreen` env var is critical — without it
PySide6 binds to the default platform plugin and tries to open an
X display / Win32 window. The conftest fixture sets it as a
fallback (`os.environ.setdefault`) but the bash run command above
sets it explicitly so a stale environment cannot break the test.

## Run all panel tests (Qt + non-Qt)

```bash
QT_QPA_PLATFORM=offscreen python -m pytest tests/panels/ -v
```

This includes the older string-grep tests in
`test_panel_ui_logic.py` and the integration tests in
`test_panel_integration.py`. Expected: ~50 passed.

## Add a new test for a new bug

The test infrastructure has 3 layers — pick whichever matches
the bug:

### Layer 1 — Widget visibility / combo items

`tests/panels/test_ih_panel_qt.py` — instantiate the real panel,
flip combos, assert `widget.isVisibleTo(panel)` and
`combo.itemText(i)`.

```python
class TestModeSwitch:
    def test_FEM_hides_air_mode(self, ih_panel):
        ih_panel._method_combo.setCurrentText("FEM")
        assert not ih_panel._widgets["air_mode"].isVisibleTo(ih_panel)
```

The `ih_panel` fixture in `tests/panels/conftest.py` builds a
fresh `IHPanel` for each test and tears it down after. No widget
state leaks between tests.

### Layer 2 — Command roundtrip (GUI -> argparse)

Same file, `TestCommandRoundtrip` class. Build the command via
`panel.build_command("model.vol")` and parse it with the
calc-script's argparse instance. Catches "GUI sends `--foo bar`
that the parser does not accept" bugs the moment they happen.

```python
def test_FEM_command_parses(self, ih_panel):
    ih_panel._method_combo.setCurrentText("FEM")
    ih_panel._widgets["workpiece_mode"].setCurrentText("SIBC")
    cmd = ih_panel.build_command("model.vol")
    parser = _calc_fem_kelvin_argparse()
    ns = parser.parse_args(cmd[2:])
    assert ns.material == "custom"
```

NOTE: The `_calc_fem_kelvin_argparse()` helper is a duplicate of
the script's argparse declaration. When you add a new CLI flag,
update BOTH the script AND the helper in the same commit.
(The duplication is intentional — we don't want the test to import
NGSolve / Radia just to validate command construction.)

### Layer 3 — Save / restore + result handler

`tests/panels/test_panel_state_restore.py` covers `save_state()`
/ `restore_state()`. Use this when you change a combo's items or
add / remove a widget.

`tests/panels/test_open_gmsh_button.py` covers the
"button enable / disable based on result dict" rule. Use this when
you add a new result key from a calc_*.py script. The mock result
dict pattern:

```python
def test_my_new_key_enables(self, ih_window, fake_msh, tmp_path):
    result = {"my_new_msh_key": fake_msh}
    assert _enable_for_result(ih_window, result, tmp_path)
```

NOTE: `_enable_for_result()` in that file replicates a slice of
`AnalysisWindow._on_finished` so the test stays focused on the
enable rule. If you change the enable rule in `radia_gui_base.py`,
update `_enable_for_result()` to match.

### Layer 4 — .vol path lifecycle and Run button enable

`tests/panels/test_run_button_browse.py` covers the
"Browse... must re-inspect the .vol so Run can enable" rule and the
QSettings round-trip path used by every interactive launch.

Use this layer when you touch:
- `AnalysisWindow._browse_vol`, `_on_vol_changed`, or the QLineEdit
  signals on `_vol_edit`
- `IHWindow.__init__` (initial `_reload_vol_info` call) or
  `set_vol_labels` semantics
- `inspect_vol_labels` (mats / bnds extraction)
- The QSettings JSON schema in `_save_settings` / `_restore_settings`

The `isolated_settings` fixture in this file monkey-patches
`radia_gui_base._SETTINGS_DIR` to a tmp dir so the developer's real
`~/.radia/radia_ih.json` is neither read (test reproducibility) nor
written (no test pollution).  When you write a save/restore test,
use this fixture or you will get flaky pass/fail depending on what
the last interactive run wrote.

Reference .vol: `tests/panels/test_3d_sibc_copper.vol` — has every
label every method needs (mats: coil/kelvin/air; bnds: sibc/source/
sink/kelvin_int/kelvin_ext).  Use it as the canonical "valid .vol
for IH validation" in new tests rather than rebuilding sample
geometries from scratch.

## Per-test fixtures

```python
def test_something(self, ih_panel):     # fresh IHPanel
def test_something(self, ih_window):    # fresh IHWindow (panel + Run/Stop/GMSH bar)
def test_something(self, qapp):         # session QApplication only
```

`ih_panel` and `ih_window` both depend on `qapp` transitively.
You don't need to add `qapp` explicitly unless your test ALSO
needs the application instance.

## Add a new panel to the test layer

When you add `radia_FOO.py` with `FOOPanel`:

1. Add a fixture to `tests/panels/conftest.py`:

   ```python
   @pytest.fixture
   def foo_panel(qapp):
       from radia_FOO import FOOPanel
       p = FOOPanel()
       yield p
       p.deleteLater()
   ```

2. Create `tests/panels/test_foo_panel_qt.py` with at least:
   - `test_default_*_is_*` for every combo's first-launch default
   - `test_*_command_parses` for every (mode, sub-mode) combination
   - widget visibility tests for every mode switch handler

3. Run the new tests + all existing panel tests:

   ```bash
   QT_QPA_PLATFORM=offscreen python -m pytest tests/panels/ -v
   ```

   Both must be green before deploying.

## CI integration (TODO)

The new tests are not yet wired into `.github/workflows/build-test.yml`.
When time permits, add a step:

```yaml
- name: PySide6 panel tests
  run: |
    pip install PySide6
    QT_QPA_PLATFORM=offscreen python -m pytest tests/panels/ -v
```

The `offscreen` platform plugin ships with PySide6 itself, no
extra system packages needed on the GitHub Actions Windows runner.

## Reference

- `tests/panels/conftest.py`            QApplication + per-test fixtures
- `tests/panels/test_ih_panel_qt.py`    IH panel widget tests
- `tests/panels/test_panel_state_restore.py`
                                        save/restore behaviour
- `tests/panels/test_open_gmsh_button.py`
                                        Open GMSH enable rule
- `tests/panels/test_run_button_browse.py`
                                        Run button enable after Browse... +
                                        QSettings round-trip (kubota 2026-05-08)
- `tests/panels/test_3d_sibc_copper.vol` reference .vol with all IH labels
- MCP knowledge: `panel_gui_pitfalls(topic="panel_qt_testing")`
- Pitfalls catalogue: `panel_gui_pitfalls()`  (full doc, 14 topics)
