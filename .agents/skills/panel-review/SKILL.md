---
name: panel-review
description: Comprehensive code review of Radia Layer-3 PySide6 panels (`src/radia/radia_*.py`). Goes beyond what panel-cli-diff (static flag matching), panel-qt-test (runtime widget wiring), and panel-preview (visual rendering) catch -- finds Layer-contract violations, solver UI duplicates, spin range vs CLI choices mismatches, subprocess error-path holes, constructor-arg vs QSettings restore-order bugs, removed-API references, and Windows cp932 encoding issues. Reports as BUG / RISK / NIT / OK with line numbers. Use AFTER editing any `radia_*.py` panel and BEFORE shipping (alongside `publish-panel`).
---

# panel-review

A comprehensive **code review** pass for Radia Layer-3 panels.  The
three existing panel skills cover specific surfaces:

- `panel-cli-diff`  — static argv ↔ argparse flag matching
- `panel-qt-test`   — runtime widget wiring (QSpinBox value reaches build_command, etc.)
- `panel-preview`   — visual rendering on real desktop Qt

This skill closes the **review gap** between them: classes of bug that
all three pass cleanly but that still ship broken behavior to users.

Real bugs caught by this skill on 2026-05-12 that the other three
missed:

| Bug class                              | Why other skills missed it |
|----------------------------------------|----------------------------|
| Duplicate solver UI entry mapping to same backend ("AMS" and "shifted AMS" both → `ams`) | cli-diff sees one flag, qt-test sees one combo, preview sees text |
| Spin range 1-3 emitted to CLI with `choices=[1,2]` | cli-diff matches flag presence, not value-range |
| `_PEEC_SOLVER_MAP` defined but never referenced | qt-test sees widget changes, cli-diff sees flag presence; nobody checks that the map is actually *read* |
| Constructor `vol_path` arg overwritten by `_restore_settings()` | qt-test exercises one path at a time |
| Module-level `from PySide6 import` with no fallback when `radia[gui]` extra missing | offscreen test imports happily; missing-extra case never tested |
| `--peec-n-peri` widget present but only forwarded by 2 of 3 builders | cli-diff DID catch this one (Silent-Drop), but only because the widget key matched |

## When to run

- AFTER any non-trivial edit to `src/radia/radia_*.py`
- AS PART OF the `publish-panel` release checklist
- BEFORE cutting a release (`release-triple` Phase 0)
- WHEN a user reports a panel bug that "should not exist" (silent-default,
  silent UI selection ignored, panel-launches-then-dies)

## What it covers (10 checks)

For each `src/radia/radia_*.py` (default: all of them — `radia_ih`,
`radia_em`, `radia_pcb`, `radia_heat`, `radia_accel`) run these
checks.  Filter to one panel with `panel-review <panel>` if needed.

### 1. Layer-3 contract violations

`radia_*.py` is **Layer 3** per AGENTS.md § Cubit Panel Architecture.
It MUST NOT import:

- `cubit` / `PyMcubit` (Layer 2 only — different Python version)
- `radia` / `ngsolve` at module top (those go through subprocess)

Grep for these imports.  Each hit is a **BUG**.

**Forbidden file formats**: Layer 3 must NEVER read, write, parse, or
expose Cubit-native artifacts in its UI or build_command paths:

- `.jou`  — Cubit journal.  Handled by Layer 1 C++ `ensure_jou_path()`
  in `RadiaComp.cpp`; the panel receives a derived `.vol` basename via
  the launcher, never the `.jou` path itself.
- `.cub5` / `.cub` — Cubit native session files.  Layer 4 reads `.vol`
  only (per AGENTS.md "Cubit/NGSolve Complete Separation Policy").

**Coil input format is solver-determined** (NOT user-selectable per
file dialog):

| Coil solver | Coil input | Widget filter MUST be |
|-------------|------------|-----------------------|
| PEEC        | `.step` / `.stp` (always — filaments derived from CAD) | `*.step *.stp` |
| BEM-A       | `.vol`     (pre-meshed surface) | `*.vol` |

This is a **hard invariant**.  PEEC mode reading a `.vol` coil, or
BEM-A reading a STEP coil, would force the wrong topology extractor
to run and is the root cause of "PEEC filament endpoint cap" and
"BEM-A source/sink load" classes of bug (see recent release notes).
Check by grepping for each coil-browse widget's `filter_str=` and
matching it against the widget's role:

```python
# Correct
add_browse("peec_step", ..., filter_str="STEP (*.step *.stp);;...")
add_browse("coil_vol",  ..., filter_str="Coil mesh (*.vol);;...")

# BUG: PEEC widget accepting .vol
add_browse("peec_step", ..., filter_str="STEP or vol (*.step *.vol)")
# BUG: BEM-A widget accepting STEP
add_browse("coil_vol",  ..., filter_str="*.vol;*.step")
```

Workpiece is always `.vol` (no STEP path).

Check pattern:

```bash
# Any of these are BUGs in radia_*.py:
grep -n '\.jou\|jou_path\|journal_file' src/radia/radia_*.py
grep -n '\.cub5\|\.cub[^i]' src/radia/radia_*.py     # exclude "cubit"
grep -n 'add_browse.*\.jou\|filter_str.*[Jj]ournal' src/radia/radia_*.py
grep -n 'QFileDialog.*\.jou' src/radia/radia_*.py
grep -n 'subprocess.*\bcubit\b' src/radia/radia_*.py
```

**Comments referencing `.jou` are OK** — they often document the
architectural boundary (e.g. "Layer 3 panels never accept .jou").  The
check is for **executable code paths** that open / save / parse the
file: `QFileDialog` filters, `add_browse(filter_str=...)`, `open(...)`
calls, basename derivation from `.jou`, etc.

**File-naming derivation**: `msh_output(vol_path, suffix)` and
`json_output(vol_path, suffix)` MUST derive their output basename from
the `.vol` path passed in, NOT from any inferred `.jou` path.  If you
see `vol_path.replace(".vol", ".jou")` or similar in the panel, that's
a BUG — the Layer 1 contract is "`.jou` basename determines all
output filenames", which is why Layer 1 saves `.jou` first and Layer 3
only sees the resulting `.vol`.

### 2. Solver-MAP / addItems consistency

A common pattern in Layer-3 panels is a `_XXX_SOLVER_MAP` dict that
translates a `QComboBox` text entry into a CLI flag value.  Two bug
classes:

**(a) UI duplicate that maps to the same backend**

```python
solver.addItems(["AMS (iterative, p=1)",
                  "shifted AMS (iterative, p=1)"])
_FEM_SOLVER_MAP = {
    "AMS (iterative, p=1)":          "ams",
    "shifted AMS (iterative, p=1)":  "shifted_ams",   # but calc_*.py rejects this
}
```

If two combo entries differ in label only but the backend choice is
identical, the second one is dead weight (and confuses users who think
they're two different solvers).

**(b) Map defined but never read**

```python
_PEEC_SOLVER_MAP = {"Dense LU (small)": ..., "HACApK (large)": ...}
# ... no reference to _PEEC_SOLVER_MAP anywhere else in the file
```

If `grep -c _XXX_SOLVER_MAP <panel>` returns 1 (just the definition),
the combo silently discards every user choice.

**(c) Map values rejected by argparse `choices=[]`**

For every `<key>: <cli_value>` pair in the map, check that
`<cli_value>` is in the target `calc_*.py` argparse `choices=[]` for
the flag the map feeds.

### 3. Spin range vs CLI `choices=[]`

```python
self.add_spin("fes_order", "Basis order:", 1, 1, 3)   # range 1-3
# but calc_inductance.py:
parser.add_argument("--h1-order", choices=[1, 2])      # rejects 3
```

For each `add_spin` whose key is passed to a calc CLI, check that the
spin's `(lo, hi)` is a **subset** of the calc's `choices=[]` (if
declared).  Mismatch = **BUG**.

If the spin range must be method-dependent (some modes accept up to 3,
others only 1-2), the panel must call `QSpinBox.setMaximum()` in the
method-change handler.  Static fixed range that's wider than any
single CLI's choices is a bug regardless of clamping at submit time.

### 4. Constructor arg vs QSettings restore order

```python
class XWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__(..., settings_key="x")
        panel = XPanel()
        self._set_panel(panel)
        # BUG: setText runs first, then _restore_settings overwrites it
        if vol_path: panel._widgets["wp_vol"].setText(...)
        self._restore_settings()
```

When a panel is launched from Cubit's Solve menu with the *current*
`.vol` path, the constructor arg should **win** over the QSettings
restore — otherwise the user opens the panel expecting today's `.vol`
and sees last week's path.

Correct pattern: `_restore_settings()` first, then override with
explicit constructor arg if non-empty.

### 5. PySide6 import guard

```python
from PySide6.QtWidgets import ...   # bare import at module top
```

`pip install radia[cubit]` (without `[gui]` extra) does NOT install
PySide6.  Cubit's Solve menu launches Layer 3 via subprocess; without
the import guard, the user clicks the menu and sees nothing —
ModuleNotFoundError on the subprocess silently kills the launch.

Required pattern:

```python
try:
    from PySide6... import ...
except ImportError as e:
    sys.stderr.write(
        "Radia <name> panel requires PySide6 ... \n"
        "  pip install --upgrade 'radia[cubit,gui]'\n")
    sys.exit(1)
```

### 6. Subprocess error path

For each `subprocess.Popen` / `subprocess.run` call in the panel:

- Is `stderr` captured (not just stdout)?
- Is non-zero exit code handled, OR is the panel relying on
  `AnalysisWindow`'s QProcess plumbing?
- If `subprocess.Popen` with detached flag (e.g. for thermal-panel
  chain launches), does the panel print a useful "Launched X" message
  but NOT claim success on the actual computation?

### 7. Output file path conventions

Per AGENTS.md "File Placement Policy", outputs (`.png`, `.msh`,
`.vtu`, JSON) must go **next to the input `.vol`**, not at repo root
or in `C:\temp\` unless the user explicitly chose a working folder.

Cross-check `msh_output(...)` and `json_output(...)` call sites
against the working-folder UX (v4.34.0+).

### 8. Removed-API references

Grep for references to removed APIs that should not appear:

- `ObjBckg([...])` array form (legacy — only the callable form is supported)
- `from radia.gmsh_builder import ...` (GmshBuilder removed)
- `CndLoop` / `CndRecBlock` / `MatSIBC` (legacy PEEC classes removed)
- `--source-mode scattered` (removed 2026-04-24)
- `sparsesolv_ngsolve` as top-level import (now `radia.sparsesolv_ngsolve`)
- `rad.FldUnits()` (removed — Radia is always meters)
- `from radia.gmsh_post_export import` is OK; flag any standalone
  `import gmsh` outside of `if __name__ == "__main__"` smoke tests

### 9. Windows cp932 encoding

Per AGENTS.md "Windows Console Encoding (cp932)", `print()` strings
must NOT contain Unicode math symbols.  Qt widget labels CAN — those
are UTF-16 internally — but anything that may reach `print()` or
stdout JSON (`json.dump` defaults to ASCII-safe) must be ASCII.

Grep the file for Unicode math symbols (`²`, `→`, `≤`, `≥`, `Δ`, `μ`,
`Ω`, `π`, `·`, `±`) and flag any in `print()`-reachable strings.
Widget labels are OK; status-line strings inside `panel_log(...)` or
`print(...)` are NIT.

### 10. Recent-commit regression patrol

Read the last ~10 commits' subject lines.  For each commit subject
that mentions a panel change ("Working folder UX", "BEM-A coil .vol
input", "PEEC filament endpoint cap fix", "spine-skips-lead
diagnostic", ...), verify the panel side of the change actually
landed:

- Was the new widget added?
- Was the corresponding CLI flag wired into the `_build_*_command`?
- Was the visibility hooked into `_on_method_changed`?

This catches the case where the calc_*.py side of a feature shipped
but the panel side regressed (or vice versa).

## How to invoke

```
panel-review              # all radia_*.py panels
panel-review radia_ih     # just IH panel
panel-review --json       # machine-readable output
```

### Delegate to a sub-agent

For the full review, delegate to a `general-purpose` agent with a
self-contained prompt that includes:

1. Path to the panel file(s)
2. Paths to the companion `calc_*.py` files
3. The 10 checks above (paste verbatim)
4. The output format below

Keep the agent's report under 800 words; ask it to verify findings by
re-reading cited lines before reporting, not just listing intent.

After the agent returns, the main session MUST verify each BUG claim
by reading the cited line and the target `calc_*.py` argparse block.
Static-analysis agents have a non-trivial false-positive rate (the
`--output` REJECT in panel-cli-diff is the canonical case —
`calc_main` in `panels/calc_common.py:1173-1177` auto-injects
`--output` if missing, so every panel's REJECT line for `--output` is
a false positive).

## Output format

```
## radia_<name>.py

### BUGS (must fix)
- line NNN: <diagnosis>. Fix: <one-line fix>

### RISKS (likely problematic)
- line NNN: <diagnosis>

### NITS (cosmetic / future)
- line NNN: <note>

### OK
- <one-sentence summary of healthy aspects>

### Verified (false positives)
- line NNN: <claim raised by sub-agent / static tool>. Why it's OK: <reason>
```

A "Verified false positive" section is important — it documents what
was checked AND ruled out, so the next reviewer doesn't re-investigate.

## Integration with other skills

`panel-review` is the **review layer**; the other panel skills are
the **regression-prevention layer**.

Recommended order before shipping a panel mode (combines with
`publish-panel`):

1. `panel-review <panel>` — find new bugs
2. Fix the BUGs and the priority RISKs
3. `panel-cli-diff <panel>` — verify CLI-flag wiring is clean
4. `panel-qt-test <panel>` — verify widget wiring at runtime
5. `panel-preview <panel>` — verify visual rendering
6. `publish-panel` — full release checklist

If `panel-review` returns BUGs, do NOT proceed to step 3.

## False-positive corpus

Document known false positives here as they're encountered, so future
runs don't re-flag them:

- `--output` REJECT in panel-cli-diff: auto-injected by
  `calc_main()` in `src/radia/panels/calc_common.py:1173-1177`.
  Waive via `# CLI-DIFF: ignore --output -- auto-injected by calc_main`.
- `--coil-msh-output` SILENT-DEFAULT in inductance modes: the panel
  passes `--msh-output` and `calc_inductance.py::_export_coil_msh_viz`
  falls back to deriving the coil-side path from the workpiece path
  (v4.25.0+).  This is by design — see memory entry
  `project_peec_viz_restored_2026_05_08.md`.

When you find a new false positive, append to this list.
