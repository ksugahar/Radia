---
name: radia-ih-panel-quality
description: End-to-end quality gate for improving and shipping the Radia induction-heating GUI panel (`src/radia/radia_ih.py`). Use when reviewing, fixing, refactoring, releasing, or deploying the IH PySide6 panel, its `calc_*.py` subprocess scripts, launcher wiring, samples, MCP docs, or first-install behavior. Also use as the template when adapting the same GUI quality workflow to other Radia panels.
---

# Radia IH Panel Quality

Run a full product-quality pass for the IH GUI panel, combining Radia's panel
skills with common GUI QA practice: HIG/usability review, visual regression,
first-run install checks, contract tests, and numeric golden validation.

This skill is an orchestration layer. When it applies, also use the existing
panel skills it names: `panel-review`, `panel-cli-diff`, `panel-qt-test`, and
`publish-panel`.

## Scope

Default target:

- Panel: `src/radia/radia_ih.py`
- Calc scripts: `src/radia/panels/calc_*.py` used by IH modes
- Shared UI base: `src/radia/radia_gui_base.py`
- Tests: `tests/panels/`
- Launcher: `src/cubit_plugin/RadiaComp.cpp`
- Public docs/knowledge: `docs/` and `packages/radia-mcp/src/radia_mcp/`

Do not widen scope to unrelated panels unless the user asks for cross-panel
generalization.

## Workflow

### 1. Establish the mode matrix

List every IH mode and submode before reading code. Include, at minimum:

- BEM / SIBC workpiece modes
- FEM / Kelvin modes
- PEEC / BEM-A / inductance-related modes if present
- Solver choices, material choices, source/sink label requirements
- Required input file kind for each mode: `.vol`, `.step`, or both

Use the matrix as the checklist for every later step. A panel mode is not
reviewed if it was not included in the matrix.

### 2. Run Radia-specific panel review

Use `panel-review radia_ih` first. Treat BUG findings as release blockers.

Verify these Radia invariants manually when the skill report is not enough:

- Layer 3 must not import Cubit, Radia physics, or NGSolve at module top.
- GUI must call headless `calc_*.py`; it must not reimplement computation.
- Constructor `.vol` argument must win over stale QSettings restore.
- PEEC STEP and BEM-A `.vol` input filters must not be mixed.
- No silent fallback, fuzzy label selection, or hidden default material path.
- Non-zero subprocess exit must surface stderr in the panel.
- Output files derive from the selected `.vol` or explicit working folder.

### 3. Check GUI-to-CLI contracts

Run:

```powershell
python tests/panels/check_panel_cli.py --panel radia_ih.py
```

If the check supports strict mode in the current repo, also run:

```powershell
python tests/panels/check_panel_cli.py --panel radia_ih.py --strict
```

Every `REJECT`, unwaived silent drop, and unwaived silent default is a bug.
Only waive a default when the panel source explains why the CLI default is the
intended user contract.

### 4. Run headless Qt behavior tests

Run the focused IH and shared panel tests:

```powershell
$env:QT_QPA_PLATFORM = "offscreen"
python -m pytest tests/panels/test_ih_panel_qt.py `
  tests/panels/test_panel_state_restore.py `
  tests/panels/test_open_gmsh_button.py `
  tests/panels/test_run_button_browse.py -v
```

For a larger pass:

```powershell
$env:QT_QPA_PLATFORM = "offscreen"
python -m pytest tests/panels/ -v
```

Add tests for any bug fixed in this pass. Prefer command-roundtrip tests that
parse the generated argv with the target calc parser surface, not string-only
assertions.

### 5. Apply general GUI QA

Perform a human-oriented review in addition to code tests:

- First launch: method combo has a valid default and no empty sections.
- Mode switching: hidden widgets do not affect commands.
- Tab order: keyboard traversal follows visual order.
- Focus: Browse/Run/error paths leave focus on the next useful control.
- Disabled states: unavailable actions explain themselves by tooltip or state.
- Error text: messages name the failed input and available labels/choices.
- DPI/font scaling: text is not clipped at 100%, 125%, and Japanese Windows fonts.
- Long paths: `.vol`, `.step`, and output paths fit or elide cleanly.
- Accessibility basics: labels are adjacent to controls and units are explicit.

Do not add in-panel explanatory prose. The panel should be operational, not a
manual embedded in a form.

### 6. Capture real desktop visual evidence

Offscreen Qt is not enough. Use `publish-panel`'s real desktop screenshot
procedure for every IH method in the mode matrix. Inspect screenshots for:

- Clipped section headers
- Overlapping rows
- Combo text truncation
- Run/Stop/Open GMSH visibility
- Orphan section headers
- Controls that jump layout when values change

If a change is visual, keep or update a reproducible screenshot command in the
report so the next reviewer can repeat it.

### 7. Validate samples and numeric golden ranges

For every user-reachable IH mode:

- Ensure a sample exists at the correct tier (`tests`, `examples`, or
  `src/radia/panels/samples`) according to AGENTS.md.
- Run the sample through the same path a user would use.
- Assert a numeric golden range, not only `status: "ok"`.
- Confirm the result JSON schema contains the keys the panel uses to enable
  Open GMSH or postprocessing buttons.

If a test is slow or validation-class, keep it outside ordinary `tests/` and
document how to run it from `examples/` or `validation/`.

### 8. First-install and deployment checks

Before calling the IH panel shipped:

- Fresh environment can import/launch the panel with `radia[cubit,gui]`.
- Missing PySide6 produces a clear install message, not a silent subprocess die.
- `cubit-plugin-install` installs the launcher and panel assets.
- Cubit launcher reaches IH on the target machine.
- The 100-machine or mdx deployment is wheel-first, not SMB file-copy surgery.
- The deployed package version matches the release or wheel under test.

Use `pyside6-health`, `radia-plugin-check`, or `deploy` when the request
includes actual machine deployment.

### 9. MCP and docs consistency

When panel behavior or user-facing options change, update public knowledge in
the same change:

- CLI flags and solver choices
- Required labels and file inputs
- Sample names and expected outputs
- Caveats and mode limitations
- Aliases users are likely to query

For Radia MCP content, update the public monorepo package under
`packages/radia-mcp/src/radia_mcp/`; do not revive retired private MCP trees.

## Output format

Use this concise report shape:

```text
IH panel quality report

Mode matrix:
- <mode>: input=<...>, calc=<...>, solver choices=<...>

Findings:
- BUG <file:line>: <issue>. Fix: <fix>
- RISK <file:line>: <issue>

Verification:
- panel-review: PASS/FAIL
- panel-cli-diff: PASS/FAIL
- panel-qt-test: PASS/FAIL
- real desktop screenshots: PASS/FAIL
- sample golden ranges: PASS/FAIL
- deploy/first-install: PASS/FAIL or not run

Docs/MCP:
- updated / already consistent / not checked
```

If no issues are found, say that clearly and name the remaining untested risk.

## Adapting to other panels

To generalize, keep the same workflow and replace only the target matrix:

- `radia_em.py`: electromagnet modes, Kelvin/periodic labels, EM samples.
- `radia_accel.py`: accelerator magnet modes and field-quality outputs.
- `radia_pcb.py`: PCB/PEEC inputs, layer/material mapping, current ports.
- `radia_heat.py`: heat-specific materials, boundary conditions, and result keys.

Do not clone this skill for every panel unless their workflow diverges. Prefer
updating this skill into a generic `radia-panel-quality` skill after IH proves
the checklist on real fixes.
