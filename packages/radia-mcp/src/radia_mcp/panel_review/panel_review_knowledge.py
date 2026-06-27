"""Knowledge topics for the Radia GUI panel review skill chain.

Each topic returns a self-contained string an AI assistant can paste
into its working context (or to a teammate) to immediately know which
check / skill to run for a class of panel bug.

Topics are versioned by content, not by version number — when a topic
needs to be updated, edit the const string here directly.  Bumping
radia-mcp's release will redistribute the latest text via PyPI.
"""

# ============================================================
# overview
# ============================================================

PANEL_REVIEW_OVERVIEW = """
# Radia GUI panel review — overview

Radia ships Layer-3 panels (`src/radia/radia_*.py`, PySide6) that
dispatch Layer-4 calc scripts (`src/radia/panels/calc_*.py`) via
subprocess.  This surface is where the most common user-visible bugs
hide: a widget that doesn't propagate to the CLI, a CLI flag that no
widget exposes, a checkbox that's always-on regardless of state, a
panel default that triggers a calc-side limitation.

The 5-skill chain in this repository covers the panel ↔ calc boundary:

| Skill | Layer | Catches |
|---|---|---|
| `panel-cli-diff` | static argv ↔ argparse | flag presence + map VALUE rejects |
| `panel-review`   | static code review (13 checks) | val()-truthy trap, removed API, Stage-1 violations, etc. |
| `panel-qt-test`  | runtime widget wiring | combo items, visibility, restore-order |
| `panel-preview`  | visual render | layout, font, clipping |
| `panel-smoke`    | END-TO-END dynamic | calc-side runtime errors triggered by panel-built command |

**When in doubt about which skill to run, use this order**:

1. `panel-review` — get a 13-check audit, fix BUGs first.
2. `panel-cli-diff` (now includes MAP-VALUE REJECT auto-check) — fix
   any REJECTs.
3. `panel-qt-test` — runtime widget regression sweep.
4. `panel-smoke` (NEW 2026-05-24) — actually run the panel and compare
   JSON against golden bands.  This is what would have caught today's
   `val()`-truthy bug and the per-element + h1-order=2 IndexError.
5. `panel-preview` — visual confirmation.
6. `publish-panel` — release checklist.

If `panel-review` returns BUGs, do NOT proceed to step 3 yet.

Topic glossary (call panel_review(topic=...)):
- `overview`            this text
- `5_skills_chain`      detailed when-to-use of each skill
- `13_checks`           full checklist used by panel-review
- `bug_catalogue`       known bug patterns, each linked to its skill check
- `val_checkbox_trap`   deep dive on the ModePanel.val() string-return bug
- `map_value_reject`    deep dive on _XXX_MAP × argparse choices mismatch
- `widget_calc_gap`     deep dive on Stage-1 (solver-switch pinning) policy
- `smoke_scenarios`     how panel-smoke runs, scenario format
- `red_flags`           quick triage table for user-reported panel issues
- `workflow`            recommended skill-ordering before shipping a panel
"""


# ============================================================
# 5_skills_chain (detailed when-to-use)
# ============================================================

PANEL_REVIEW_5_SKILLS_CHAIN = """
# 5-skill chain — when to use each

## panel-cli-diff (static)
**File**: `tests/panels/check_panel_cli.py`
**Invoke**: `python tests/panels/check_panel_cli.py [--panel radia_ih.py] [--strict]`
**Output**: per-method table of OK / REJECT / SILENT-DEF / MAP-VALUE REJECT / ORPHAN flags.

Catches:
1. Panel emits `--foo` but no calc declares it → REJECT
2. Calc has `--bar` but no panel emits → SILENT-DEF
3. Panel defines widget X that nothing reads → ORPHAN
4. (NEW 2026-05-24) `_XXX_MAP` leaf value not in any calc argparse
   `choices=[]` → MAP-VALUE REJECT

Use AFTER editing any `radia_*.py` panel OR `calc_*.py` argparse block.

## panel-review (static, 13 checks)
**File**: `.claude/skills/panel-review/SKILL.md`
**Invoke**: launch via Skill tool, or delegate to general-purpose agent
with the 13 checks pasted verbatim.

Catches:
1. Layer-3 contract violations (top-level cubit/ngsolve import, .jou
   read in executable code path, coil-input file format vs solver
   mismatch)
2. Solver-MAP / addItems consistency (UI duplicate, dead map, choices
   mismatch)
3. Spin range vs CLI choices
4. Constructor arg vs QSettings restore order
5. PySide6 import guard
6. Subprocess error path
7. Output file path conventions
8. Removed-API references (FldUnits, ObjBckg array, gmsh_builder, etc.)
9. Windows cp932 encoding (Unicode in print())
10. Recent-commit regression patrol (calc shipped, panel forgot)
11. (NEW) val() checkbox truthy-trap (`if self.val("ck"):` always True)
12. (NEW) solver-map × argparse choices automated cross-check
13. (NEW) widget default vs known calc-side limitation (e.g. P2 + per-DOF)

Use AFTER any non-trivial edit to `radia_*.py`.

## panel-qt-test (runtime widget)
**File**: `tests/panels/test_*_qt.py` + `panel-qt-test` skill
**Invoke**: `pytest tests/panels/test_*_qt.py`

Catches:
- combo items match expected list
- visibility rules per method
- save/restore round-trip
- spin range clamp behaviour

Use AFTER widget visibility or save/restore changes.

## panel-preview (visual)
**File**: `.claude/skills/panel-preview/preview.py`
**Invoke**: `python .claude/skills/panel-preview/preview.py --panel radia_ih`

Catches:
- widget actually appears in the right method
- layout clipping
- ordering of widgets within sections

Use when "show me what the panel looks like" or to verify a new widget
visually appeared.

## panel-smoke (END-TO-END dynamic, NEW 2026-05-24)
**File**: `.claude/skills/panel-smoke/smoke.py` + `tests/panels/smoke_scenarios.py`
**Invoke**: `python .claude/skills/panel-smoke/smoke.py [--id <scenario>]`

Catches **everything the others cannot** because it runs the actual
subprocess that the panel build_command produces, then validates the
JSON output against a golden band:

- `val()`-truthy: panel-built command would emit `--esim-per-panel`
  even when checkbox is unchecked → scenario expects `False`, gets
  `True` → mismatch FAIL
- P2 + per-DOF IndexError: subprocess JSON has `{"error": "..."}`
  field → scenario FAIL
- Numerical drift: any calc-layer change that shifts P_wp > tol → FAIL
- Panel default triggers calc limitation → FAIL

Scenarios are catalogued in `tests/panels/smoke_scenarios.py`.  Reference
scenarios that MUST stay green:
- `ih_peec_bem_scalar_I100`           — P_wp ≈ 30.51 W (sweep_v2 scalar)
- `ih_peec_bem_per_element_I100`      — P_wp ≈ 18.75 W (sweep_v2 per-panel, IGTE Fig. 1)

Add a scenario for every new panel mode you ship.
"""


# ============================================================
# 13_checks (full checklist)
# ============================================================

PANEL_REVIEW_13_CHECKS = """
# panel-review skill — 13 checks (2026-05-24 + 3 additions)

For each `src/radia/radia_*.py` panel, run these checks.

## 1. Layer-3 contract violations
- No top-level `import cubit` / `import radia` / `import ngsolve`
- No `.jou` / `.cub5` read in executable code paths
- Coil input format matches solver:
    - PEEC → STEP only (`*.step *.stp`)
    - BEM-A → .vol only (`*.vol`)
- Workpiece always `.vol` (no STEP path)

## 2. Solver-MAP / addItems consistency
- UI duplicate that maps to same backend (dead label)
- Map defined but never read in any build_command
- Map values rejected by target argparse `choices=[]`

## 3. Spin range vs CLI choices
For each `add_spin(key, ..., lo, hi)` whose key feeds a CLI flag with
declared `choices=[]`, verify `(lo, hi)` is a subset of choices.

## 4. Constructor arg vs QSettings restore order
`_restore_settings()` MUST run before constructor-arg override of any
widget (so today's vol_path wins over last week's saved state).

## 5. PySide6 import guard
Bare `from PySide6 import` at module top crashes silently when user
installs `radia[cubit]` without `[gui]`.  Wrap in try/except + stderr
hint.

## 6. Subprocess error path
- stderr captured (not just stdout)
- Non-zero exit handled OR delegated to AnalysisWindow base
- Detached Popen prints "Launched X" but does NOT claim success

## 7. Output file path conventions
`msh_output(vol_path, suffix)` and `json_output(...)` derive from the
`.vol` path passed in, never from any `.jou` path.

## 8. Removed-API references
Grep for:
- `ObjBckg([...])` array form (legacy)
- `gmsh_builder` (removed)
- `CndLoop` / `MatSIBC` (legacy PEEC removed)
- `--source-mode scattered` (removed 2026-04-24)
- `sparsesolv_ngsolve` top-level (now `radia.sparsesolv_ngsolve`)
- `rad.FldUnits()` (removed)
- top-level `import gmsh` outside `if __name__ == "__main__"`

## 9. Windows cp932 encoding
Unicode math symbols (², →, ≤, Δ, μ, Ω, π, ·, ±) MUST NOT appear in
`print()` strings.  Qt widget labels are fine (UTF-16 internally).

## 10. Recent-commit regression patrol
Read `git log --oneline -10 <panel>`.  For each commit subject that
mentions a panel-side feature, verify the panel actually got the
widget + emit.

## 11. val() checkbox truthy-trap (NEW)
`ModePanel.val(key)` returns the STRING `"1"` / `"0"` for QCheckBox
(both truthy).  `if self.val("ck"):` is always True.

Correct: `if self._widgets["ck"].isChecked():`

## 12. Solver-map × argparse choices cross-check (NEW)
Now automated in panel-cli-diff Check 4 (MAP-VALUE REJECT).

## 13. Widget default vs known calc-side limitation (NEW)
Known per-DOF ESIM constraint (v4.67-v4.72):
- `extract_H_t_per_dof_grad` expects P1 vertex DOFs.  Panel must
  clamp `fes_order` to 1 when `esim_per_panel` is on.

For every newly added solver-switch widget, audit visibility AND
default vs method-specific calc-side constraints.
"""


# ============================================================
# bug_catalogue (the patterns from today)
# ============================================================

PANEL_REVIEW_BUG_CATALOGUE = """
# Known panel bug patterns + detection skills

## Pattern A: val() checkbox truthy-trap
**Symptom**: A CLI flag (`--esim-per-panel`, `--newton`, `--require-X`)
is emitted regardless of the corresponding checkbox state.
**Root cause**: `self.val("checkbox_key")` returns `"1"` or `"0"`
(both truthy).  `if self.val(...):` is always True.
**Detection**: panel-review Check 11, panel-smoke (scenario expects
checkbox=False but JSON output has the feature ON).
**Fix**: `if self._widgets["key"].isChecked():`  (see
`_heat_panel.py:420` for the reference idiom)

## Pattern B: AMG vs AMS map-value typo
**Symptom**: User selects "AMG (Compact)" → subprocess exits 2 with
`error: argument --solver: invalid choice: 'amg'`.
**Root cause**: `_SOLVER_MAP["AMG (Compact)"] = "amg"` but
calc's argparse has `choices=["...","ams"]` — no "amg" entry.
**Detection**: panel-cli-diff Check 4 (MAP-VALUE REJECT).
**Fix**: rename the UI label to match the backend OR remove the dead
choice from the map.

## Pattern C: Widget-vs-calc gap (Stage-1 violation)
**Symptom**: User can't reproduce a paper-cited / release-noted feature
through the panel — has to drop to CLI.
**Root cause**: calc layer shipped the flag in 1 commit, panel layer
was never updated.  Per CLAUDE.md "Panel Design Workflow Policy"
Stage 1 the solver-switch MUST be pinned in the panel.
**Detection**: panel-cli-diff SILENT-DEF row + panel-review Check 10
(recent-commit regression patrol), panel-smoke (paper scenario fails).
**Fix**: add widget(s) + wire build_command emission + visibility.
For checkbox: use `isChecked()` not `val()` (see Pattern A).

## Pattern D: Widget default × calc-side limitation
**Symptom**: Panel-built command exits 0 but JSON has
`{"error": "index ... out of bounds"}` or similar internal calc fault.
**Root cause**: Panel default (e.g. `fes_order=2`) hits a calc-layer
limitation (e.g. per-DOF extractor is P1-only) without the panel
auto-clamping.
**Detection**: panel-smoke (JSON has `error` key).
**Fix**: clamp the spinbox max in the relevant `_on_X_changed`
handler (see `radia_ih.py:835-844` for the existing
method-dependent pattern, applied at `_on_impedance_changed` for the
per-DOF case).

## Pattern E: stateChanged → handler-cascade recursion
**Symptom**: Panel init throws RecursionError because handler A
fires handler B which fires handler A.
**Root cause**: connecting a checkbox's stateChanged to a handler
that calls a method handler that re-fires impedance handler.
**Detection**: panel-qt-test fixture instantiation raises.
**Fix**: replicate the needed logic in-place rather than recursing
back into another handler.

## Pattern F: Stale test assertions
**Symptom**: Two pre-existing test failures appear in any panel test
run, but no code regression — assertion text matches an older panel
state.
**Detection**: pytest -v shows assertion mismatch on a *count* or
*list* of widgets.
**Fix**: update the assertion to match current state, add a comment
noting which feature change made the old assertion stale.

## Pattern F2: Legacy .vol row tests after working-folder migration
**Symptom**: IH panel tests fail with missing
`AnalysisWindow._on_vol_changed`, `_browse_vol`, or stale `_vol_edit`
assertions.
**Root cause**: radia 4.35.0 moved file selection to a top-level
`working_folder` plus panel-owned browse rows.  IH now owns its
workpiece mesh at `IHPanel._widgets["wp_vol"]`, and label inspection
runs through `IHPanel._on_wp_vol_changed_text`.
**Detection**: `tests/panels/test_run_button_browse.py` should exercise
`panel.wp_vol_path()`, saved `working_folder` + relative `panel.wp_vol`,
and constructor `vol_path` override semantics.
**Fix**: update tests and review checklists to the panel-owned
`wp_vol` contract.  Do not resurrect the hidden `_vol_edit` shim as a
source of truth; it exists only for old settings JSON compatibility.

## Pattern G: Layout regression in offscreen render
**Symptom**: panel-preview screenshot shows label overlap on a
section header.
**Root cause**: dynamic status label placed at the row boundary of a
section header.
**Caveat**: offscreen Qt renderer compresses fonts, so PNG visual
overlap might be a rendering artefact, NOT a real bug.  Check
real-Qt render (panel-preview default) or .txt widget list before
filing.

## Pattern H: Other-user lock blocking PDF / file overwrite
**Symptom**: `cp` to W:\\ fails with "Device or resource busy" even
after killing all known PDF holders.
**Root cause**: file is open in another user's RDP session on
100号機 (shared W:\\ share).
**Fastest fix**: ask the file owner (e.g. kubota) to sign out,
rather than spelunking through cross-session process kill.
"""


# ============================================================
# val_checkbox_trap deep dive
# ============================================================

PANEL_REVIEW_VAL_CHECKBOX_TRAP = """
# val() checkbox truthy trap — deep dive

`src/radia/radia_gui_base.py:474-475`:

```python
elif isinstance(w, QCheckBox):
    return "1" if w.isChecked() else "0"
```

Both `"1"` and `"0"` are TRUTHY in Python.  So:

```python
# BUG (silent): emits --foo unconditionally
if self.val("enable_foo"):
    cmd.append("--foo")
```

is equivalent to `if True:` because `bool("0") == True`.

Correct idiom (matches `_heat_panel.py:420`):

```python
if self._widgets["enable_foo"].isChecked():
    cmd.append("--foo")
```

Alternative explicit-string idiom (works but more fragile if the
return convention changes):

```python
if self.val("enable_foo") == "1":
    cmd.append("--foo")
```

## How to detect

panel-review Check 11 grep:

```
grep -nE 'if self\\.val\\("([a-z_]+)"\\):' src/radia/radia_*.py
```

For each match, cross-reference the key against
`self.add_check("<key>", ...)` sites in the same file.  If the key is
a QCheckBox AND the call is NOT inside an explicit `== "1"` comparison,
flag as RISK.

panel-smoke catches it dynamically: scenario sets checkbox=False but
JSON output has the feature flag on.

## Why we ship the string convention at all

`val()` returns are usually passed directly into the subprocess argv
list as string args (`cmd += [..., self.val("freq"), ...]`), so
stringifying at the source kept the call sites uniform.  Checkbox is
the one widget type where the user almost always wants a boolean for
control-flow, not a string for CLI emission.  Until `val()` is
refactored to return bool for QCheckBox, use `isChecked()` directly
for any control-flow site.
"""


# ============================================================
# map_value_reject deep dive
# ============================================================

PANEL_REVIEW_MAP_VALUE_REJECT = """
# Map-value REJECT — deep dive

## The bug class

A panel commonly defines:

```python
_SOLVER_MAP = {
    "Direct (PARDISO)": "pardiso",
    "AMG (Compact)":    "amg",       # <-- BUG: "amg" is not in calc choices
    "AMS (Compact)":    "ams",
    "BDDC":             "bddc",
}
```

If `calc_*.py` declares:

```python
parser.add_argument("--solver", choices=["auto", "pardiso", "bddc",
                                          "iccg", "ams"])
```

then selecting "AMG (Compact)" in the panel sends `--solver amg` →
subprocess exits 2 with `error: argument --solver: invalid choice:
'amg' (choose from auto, pardiso, bddc, iccg, ams)`.

The pre-2026-05-24 panel-cli-diff only checked flag *presence*, not
*values*, so this slipped through.

## Detection (automated, panel-cli-diff Check 4)

`tests/panels/check_panel_cli.py` Check 4 (MAP-VALUE REJECT) scans
every panel for class-level `_XXX_MAP` dicts, extracts leaf string
values, and checks each against the UNION of `choices=[...]` literals
across all `calc_*.py` argparse calls.  Values that are CLI-flag-
shaped (lowercase, no spaces, no parens) but absent from any
choices list are flagged as REJECT.

Sample output (pre-fix):

```
=== src/radia/radia_em.py ===
  MAP-VALUE REJECT ( 1): _SOLVER_MAP['amg']
```

## Waiver mechanism

For map values that are consumed internally (not as `--flag <value>`
CLI args), add to the panel source:

```python
# CLI-DIFF: ignore-map-value <value1> <value2> -- reason
_XXX_MAP = {
    "Label A": "internal-state-value",  # not a CLI flag value
}
```

The checker will skip the listed values.

## Fix template

Rename UI label to match the backend label, then update the map:

```python
# Before
_OMEGA_SOLVERS = ["Direct (PARDISO)", "AMG (Compact)", "BDDC"]
_SOLVER_MAP = {..., "AMG (Compact)": "amg", "AMS (Compact)": "ams", ...}

# After (today's fix in radia_em.py)
_OMEGA_SOLVERS = ["Direct (PARDISO)", "AMS (Compact)", "BDDC"]
_SOLVER_MAP = {..., "AMS (Compact)": "ams", ...}  # AMG entry removed
```

## Related

panel-review Check 12 documents this class and refers to the
automated check.
"""


# ============================================================
# widget_calc_gap deep dive
# ============================================================

PANEL_REVIEW_WIDGET_CALC_GAP = """
# Widget-vs-calc gap (Stage-1 violation) — deep dive

## The CLAUDE.md policy this enforces

> Stage 1 — Enumerate the app-specific variables.  Write down every
> knob the user of this specific application might want to change.
> ...  Pin the solver-specific variables (mesh size, frequency,
> material, source current, ...) and the **solver-switch variable**
> itself (e.g. `--impedance-model linear|esim`,
> `--solver pardiso|ams`).

When a developer adds a CLI flag to a calc script (Stage 2 = Layer 4)
without surfacing it in the panel (Stage 3 = Layer 3), users cannot
reproduce the feature through the UI even though the calc is fully
capable.  For paper-headline contributions this is especially bad —
the paper's whole value depends on users being able to reproduce.

## Real example (IGTE 2026 paper)

| commit | calc-side change | panel-side change |
|---|---|---|
| `630527d4` | `fix(esim): per-DOF |H_t| via triangle-gradient` | none |
| `c383cd5e` | `esim(anderson): safeguarded Anderson` | none |
| `0c401f03` | `esim(anderson): Anderson-type-II acceleration` | none |
| `36ef13fc` | `feat(calc_inductance): emit esim_per_panel_H_t array` | none |

Three calc commits in May 2026 added `--esim-per-panel` +
`--esim-anderson-m`, but `radia_ih.py` was never touched.  Result:
the IGTE paper's headline (per-element ESIM gives -38.5% disagreement
vs scalar) was UNREACHABLE from the panel UI for ~6 weeks before
kubota noticed.

## How to detect

(a) **panel-cli-diff SILENT-DEF row**: every flag in the calc that no
panel emits shows up.  But SILENT-DEF can be intentional (e.g. an
expert-mode flag); requires human judgment to triage.

(b) **panel-review Check 10** (recent-commit regression patrol): read
`git log --oneline -10 <calc>` and `git log --oneline -10 <panel>`.
If calc had feature commits but panel had none, flag.

(c) **panel-smoke scenario for the paper-headline case**: if the
scenario can't be expressed in the panel widgets, the gap is
self-evident.

## Fix template

The IGTE 2026 patch (2026-05-24, radia_ih.py) is the reference:

```python
# In _build_ui, after existing ESIM widgets:
per_panel_check = self.add_check(
    "esim_per_panel",
    "Per-DOF ESIM (resolves saturation hot-spots; IGTE 2026 paper)",
    default=False)
per_panel_check.stateChanged.connect(
    lambda _: self._on_impedance_changed(self.val("impedance_model")))
self.add_spin("esim_anderson_m", "Anderson memory m:", 5, 0, 20)
self.add_line("esim_relax", "Karl relax alpha:", "0.5")

# In _on_impedance_changed visibility loop:
for key in ("bh_file", "esim_max_iter", "esim_tol",
            "esim_per_panel", "esim_anderson_m", "esim_relax"):
    self._set_row_visible(key, is_esim)

# In _build_X_command (where X has ESIM support):
if self._widgets["esim_per_panel"].isChecked():    # <-- isChecked() NOT val()
    cmd.append("--esim-per-panel")
cmd += ["--esim-anderson-m", str(self.val("esim_anderson_m")),
        "--esim-relax",       self.val("esim_relax")]
```

NB: not every method accepts every flag.  Verify per-method calc-side
support and hide widgets / skip emission accordingly (e.g.
`calc_fem_kelvin.py` accepts `--esim-per-panel` but NOT
`--esim-anderson-m`).
"""


# ============================================================
# smoke_scenarios deep dive
# ============================================================

PANEL_REVIEW_SMOKE_SCENARIOS = """
# panel-smoke scenarios — format + workflow

## Why scenarios

panel-smoke (NEW 2026-05-24) drives a panel via offscreen Qt with a
specific dict of widget values, calls `build_command`, runs subprocess,
and validates the JSON against a golden band.  This is the only skill
that catches:

- `val()` checkbox truthy bug (scenario expects checkbox=False but
  JSON has feature ON)
- panel default that triggers calc IndexError (JSON has error key)
- numerical drift > tolerance (regression on the headline)

## Scenario file

`tests/panels/smoke_scenarios.py` — list of dicts:

```python
{
    "id": "ih_peec_bem_per_element_I100",
    "panel": "radia_ih",
    "window_class": "IHWindow",
    "method_const": "METHOD_PEEC_BEM",
    "vol": str(SAMPLES / "ih_bem_sample_p1.vol"),
    "set": {
        "impedance_model": "Nonlinear ESIM (experimental, WIP)",
        "peec_step": str(SAMPLES / "ih_fem_kelvin_demo_coil.step"),
        "bh_file":   str(SAMPLES / "em_sample_bh.txt"),
        "freq":            "50000",
        "current":         "100.0",
        "wp_sigma":        "2e6",
        "mu_r":            "100",
        "half_thickness":  "0.005",
        "esim_per_panel":  True,        # bool for QCheckBox
        "esim_anderson_m": 5,            # int for QSpinBox
        "esim_max_iter":   30,
        "esim_relax":      "0.5",        # str for QLineEdit
    },
    "builder": "_build_peec_bem_command",
    "golden": {
        "P_wp_W":          {"value": 18.75, "tol_rel": 0.02},
        "L_coil_nH":       {"value": 100.67, "tol_rel": 0.01},
        "esim_iterations": {"value": 7, "tol_abs": 3},
        "esim_converged":  {"value": True},
        "esim_per_panel":  {"value": True},
        "esim_anderson_m": {"value": 5},
    },
    "timeout_s": 600,
}
```

## Invocation

```bash
python .claude/skills/panel-smoke/smoke.py                    # all
python .claude/skills/panel-smoke/smoke.py --list             # IDs
python .claude/skills/panel-smoke/smoke.py --id <id> -v       # one verbose
python .claude/skills/panel-smoke/smoke.py --panel radia_ih   # by panel
```

## Adding a new scenario

When you add a panel mode or a new sample:

1. Run the panel by hand once (panel-preview can help verify visual);
   capture the JSON your widgets produce.
2. Append a scenario dict with the values you set + a golden table
   matching the JSON keys you want to lock.
3. Use `tol_rel` for floating-point numerics that may drift with
   solver internals (typical: 1-5%).  Use `tol_abs` for iteration
   counts that round to integers.
4. Use exact-match for booleans (`{"value": True}` only).
5. Re-run `--id <new_id>` to confirm.

## Golden source of truth — IH

- `examples/ih_esim_benchmark/sweep_data/sweep_results.json` is the
  paper-citable reference for IH PEEC+BEM at all 32 (I, f) cells.
- For new IH scenarios, derive golden from sweep_v2 or run the
  scenario manually to capture a reference JSON before adding it.

## Reference scenarios (must stay green at every release)

- `ih_peec_bem_scalar_I100`      — sweep_v2 scalar (P_wp ≈ 30.51 W)
- `ih_peec_bem_per_element_I100` — sweep_v2 per-panel (P_wp ≈ 18.75 W,
                                    IGTE 2026 Fig. 1 representative)

Add scenarios for: EM Kelvin Benchmark, PCB PEEC, motor lamination
(when widget catalogue extension lands).
"""


# ============================================================
# red_flags (triage table)
# ============================================================

PANEL_REVIEW_RED_FLAGS = """
# Red-flag triage — user reports "panel does something weird"

| User report | First suspect | Skill to run |
|---|---|---|
| "I clicked Run and nothing happened" | PySide6 import guard missing OR subprocess detach silent crash | panel-review Check 5 + 6 |
| "Wrong solver is being used regardless of what I select" | val() checkbox-truthy OR _SOLVER_MAP dead | panel-review Check 11 + 2(b), panel-smoke |
| "Subprocess exits 2 with argparse error" | panel emits an unknown flag OR map value rejected | panel-cli-diff (REJECT + MAP-VALUE REJECT) |
| "Result is numerically wrong by ~factor" | calc-layer regression OR widget default mismatch | panel-smoke against golden |
| "I can't reproduce a paper / release-note result from the panel" | widget-vs-calc gap (Stage-1 violation) | panel-cli-diff SILENT-DEF + panel-review Check 10 + 13 |
| "Panel opens with old vol_path even though Cubit launched it with today's" | constructor arg overwritten by QSettings restore | panel-review Check 4 |
| "JSON output has {'error': 'index N is out of bounds for size N'}" | widget default × calc-side limitation (e.g. P2 + per-DOF) | panel-review Check 13, panel-smoke |
| "Status label overlaps a section header in the screenshot" | possibly real layout bug OR offscreen-Qt artefact | panel-preview real-Qt + .txt cross-check |
| "Can't write the PDF — Device or resource busy" | other-user RDP session holds the file on 100号機 | ask owner to sign out (not a code bug) |
| "Panel widget I changed isn't reaching the calc" | widget-key typo OR build_command forgot to emit | panel-cli-diff SILENT-DEF / ORPHAN |
| "Panel tests fail after I edited nothing" | stale assertion (test drifted from current panel state) | check git blame on the assertion, update |

## When in doubt: run panel-smoke first.

It exercises the same code path the user did and prints PASS/FAIL on
named golden fields, often pointing directly at the root cause in
~1-10 minutes of subprocess time.
"""


# ============================================================
# workflow (recommended skill ordering)
# ============================================================

PANEL_REVIEW_WORKFLOW = """
# Recommended workflow — before shipping a panel mode

## After editing a Layer-3 panel (radia_*.py)

1. **panel-review** — fix all BUGs; check RISKs.  If BUGs remain, stop
   here and fix them first.
2. **panel-cli-diff** — fix REJECT + MAP-VALUE REJECT rows.
3. **panel-qt-test** — `pytest tests/panels/test_<panel>_qt.py`
4. **panel-smoke** — `python .claude/skills/panel-smoke/smoke.py
   --panel <panel>`  (add a new scenario if your mode is not yet
   catalogued)
5. **panel-preview** — visual confirmation that widgets render in
   the right place
6. **publish-panel** — release checklist (also runs panel-preview
   + panel-qt-test as gates)

## After editing a Layer-4 calc script (calc_*.py)

1. **panel-cli-diff** — does your new flag have a panel widget?  If
   SILENT-DEF on a paper-headline flag, that's a widget gap.
2. **panel-review Check 10** — recent-commit regression patrol:
   panel got the corresponding widget?
3. **panel-smoke** — does the affected scenario still pass?

## Before a release-qud

`release-qud` Phase 0 should include:

```bash
python tests/panels/check_panel_cli.py --strict
python .claude/skills/panel-smoke/smoke.py
pytest tests/panels/test_*_qt.py -q
```

All three must pass before tagging.

## When fielding a user-reported bug

See `panel_review(topic='red_flags')` for the triage table.  Default:
run `panel-smoke --id <closest scenario> -v` first; the verbose
output points at the failing field within seconds.
"""


# ============================================================
# Dispatch table
# ============================================================

TOPICS = {
    "overview":           PANEL_REVIEW_OVERVIEW,
    "5_skills_chain":     PANEL_REVIEW_5_SKILLS_CHAIN,
    "13_checks":          PANEL_REVIEW_13_CHECKS,
    "bug_catalogue":      PANEL_REVIEW_BUG_CATALOGUE,
    "val_checkbox_trap":  PANEL_REVIEW_VAL_CHECKBOX_TRAP,
    "map_value_reject":   PANEL_REVIEW_MAP_VALUE_REJECT,
    "widget_calc_gap":    PANEL_REVIEW_WIDGET_CALC_GAP,
    "smoke_scenarios":    PANEL_REVIEW_SMOKE_SCENARIOS,
    "red_flags":          PANEL_REVIEW_RED_FLAGS,
    "workflow":           PANEL_REVIEW_WORKFLOW,
}


def get_panel_review_documentation(topic: str = "overview") -> str:
    """Return the requested topic, or a topic list if 'all' / unknown."""
    if topic == "all":
        return "\n\n".join(f"# Topic: {k}\n{v}" for k, v in TOPICS.items())
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic: {topic!r}\n\n"
        f"Available topics:\n"
        + "\n".join(f"  - {k}" for k in TOPICS)
    )
