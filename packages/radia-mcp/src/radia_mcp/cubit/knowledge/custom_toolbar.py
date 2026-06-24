"""
Custom Toolbar / In-Cubit PySide6 knowledge base.

Covers the in-process Cubit Python + PySide6 workflow for building custom
toolbars, dialogs, and encapsulated multi-step workflows.

Sources:
  - Coreform webinar "How to create a custom GUI in Coreform Cubit (Intro)"
    https://www.youtube.com/watch?v=TnZyZHDMOwA  (Carl McKelvey)
  - Coreform Cubit 2025.12 / LAB deployment (PySide6 shipped in distribution).
  - Coreform examples repo: https://github.com/coreform-llc
    (DAGMC toolbar, tire cross-section toolbar -- both .tar.gz packaged).

NOTE: This is *different* from two adjacent topics:
  (a) the Radia-NGSolve analysis panels described in `panel_conventions.py`
      -- those are standalone PySide6 apps running in Python 3.12 and
      launching Cubit as a separate process.  The knowledge here covers
      PySide6 code that runs *inside* Cubit's embedded Python 3.10
      interpreter, bound to Cubit's own Qt main window.
  (b) the C++ SDK plugin approach described in `cpp_sdk.py` -- that is for
      deep customization (new command panels with navigation nodes, new
      command-line verbs with hotkeys, deep menu surgery) where Python
      cannot reach.  Use the Python toolbar path here whenever it
      suffices; reach for the C++ SDK only when a hotkey-bound new
      command or a fully custom navigation-node command panel is needed.
"""

TOOLBAR_OVERVIEW = """
# In-Cubit Custom Toolbar Overview (Coreform Cubit 2025.12+)

Coreform Cubit 2025.12 ships PySide6 as part of its distribution.
You can extend the Cubit GUI with custom toolbars that trigger
Cubit command-language scripts, Python scripts, or nested command panels.

## When to use this

- Repetitive multi-step workflows (e.g., "clean geometry -> imprint+merge
  -> block assign -> mesh -> export") that users want as one button.
- Parameter-driven variants where a small dialog collects values before
  running a script.
- Team-shared workflows (packaged `.tar.gz`, distributed + imported via
  right-click on the toolbar panel).

## When NOT to use this

- If the GUI must run outside Cubit (e.g., a launcher that spawns Cubit as
  a subprocess, or uses Python 3.12 for NGSolve), use the Radia-NGSolve
  panel convention instead (`cubit_docs topic=panel_conventions`). Those
  run in a separate interpreter and pipe commands to Cubit via journal
  files.

## Entry points

| How to open the editor         | Path                                            |
|--------------------------------|-------------------------------------------------|
| Menu                           | `Tools -> Custom Toolbar Editor`                |
| Toolbar                        | Hammer icon on the main toolbar                 |

## Button types

| Type                   | When to pick                                     |
|------------------------|--------------------------------------------------|
| Simple tool button     | Cubit command-language script (no GUI needed)    |
| Python script button   | Python (with or without PySide6 dialog)          |
| Command panel opener   | Launches another command panel (for grouping)    |

## Reference toolbars (Coreform GitHub)

- `coreformLLC/cubit-toolbars` — DAGMC (neutronics, Monte Carlo) + tire
  cross-section workflow. Study these before writing your own.
- Packaging scripts (`package-win.ps1`, `package-linux.sh`) produce
  `.tar.gz` archives. Import via right-click -> Import in the toolbar
  panel; archive name becomes the toolbar namespace.

## License posture

- Coreform toolbar examples are MIT-licensed.
- PySide6 is LGPL. Distributing closed-source toolbars that dynamically
  link PySide6 is fine under LGPL; if you need to hide business logic,
  factor it into a C++ component (see Coreform Norbert's calculus sample)
  rather than trying to obfuscate Python.

## Relation to other modules

- `panel_conventions_knowledge.py` - external PySide6 launcher (Python
  3.12, NGSolve, spawns Cubit as subprocess).
- `scripting_knowledge.py` - Cubit command-language + Python API
  reference. The toolbar scripts described here *use* those APIs.
"""

TOOLBAR_PYTHON_SCRIPT_CONVENTIONS = """
# Cubit In-Process Python Script Conventions (CRITICAL)

Scripts executed from a custom toolbar run inside Cubit's embedded Python
interpreter (currently Python 3.10 for Cubit 2025.12 / LAB). This interpreter has
subtle differences from standard CPython that break naive scripts.

## Required script header

```python
#!python
# The shebang tells Cubit "this is a Python script" (not a journal file).
# Without it the script is sometimes parsed as Cubit command language.

# The `cubit` module is ALREADY IMPORTED by the host. Do NOT write
# `import cubit` at the top of a toolbar script - Cubit's importer can
# fail on the name shadowing. Just use `cubit.cmd(...)` directly.

# When a script is run via the toolbar, __name__ is '_coreform_cubit'
# (not '__main__'). Use this to detect in-Cubit execution:
if __name__ == '_coreform_cubit':
    # in-Cubit entry point
    dlg = MyDialog()
    dlg.exec()
```

## Import statement gotcha (HIGH IMPACT)

Cubit's Python importer **chokes on parenthesized multi-line imports**.
This is a known limitation called out by Coreform in their 2025 GUI
webinar. Use backslash continuations instead.

```python
# WRONG - Cubit's in-process importer often fails on this:
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QLineEdit,
    QLabel,
    QDialogButtonBox,
)

# CORRECT - backslash continuations:
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLineEdit, \\
    QLabel, QDialogButtonBox

# Or single-line (if it fits):
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLineEdit, QLabel, QDialogButtonBox
```

Linter rule: `cubit-toolbar-import-parens` flags parenthesized imports in
scripts that look like toolbar scripts (detected by `#!python` shebang,
`_coreform_cubit` namespace check, or location under a `scripts/` folder).

## Float conversion with error UX

LLM-generated dialogs often call `float(line_edit.text())` unguarded.
Always wrap in try/except and show a QMessageBox so the user gets
feedback instead of a silent stack trace in the Cubit terminal.

```python
try:
    size = float(self.size_edit.text())
    radius = float(self.radius_edit.text())
except ValueError as exc:
    QMessageBox.warning(self, "Invalid input",
                        f"Enter a numeric value: {exc}")
    return
cubit.cmd(f"create brick x {size}")
cubit.cmd(f"create cylinder radius {radius} height {size * 1.5}")
cubit.cmd("subtract volume 2 from volume 1")
```

## String templating vs object API

The `cubit.cmd("...")` string templating is the workhorse - it is the
most widely used and most thoroughly tested interface, because Cubit's
own GUI drives core-Cubit through the same Python binding layer (wrapped
via SWIG).

A Pythonic object API exists (e.g., `cubit.fire_array(...)` for raycast
queries), but it is less maintained and mostly read-only. Do not try to
"abstract away" the command language; use f-strings + numpy instead.

```python
import numpy as np
# Generate commands with numpy + f-strings, then issue them
for i, x in enumerate(np.linspace(0, 10, 5), start=1):
    cubit.cmd(f"create vertex {x} 0 0")
```
"""

TOOLBAR_CLARO_PARENT_HELPER = """
# Claro Parent Helper (prevents dialog-behind-window problem)

When a PySide6 dialog is created without a parent inside Cubit, it can
fall behind the main Cubit window and disappear (especially on Windows
with multiple monitors). The fix is to find Cubit's QMainWindow (named
`claro`) and use it as the dialog parent.

This is the standard helper Coreform demonstrates in the custom GUI
webinar and recommends copying into your toolbar's utility script.

## Helper function

```python
from PySide6.QtWidgets import QApplication, QMainWindow

def find_claro():
    \"\"\"Return Cubit's main window (QMainWindow named 'claro'), or None.

    Call this once at script start, then pass the result as `parent=`
    to every QDialog you create. Prevents Windows from losing the
    dialog behind the Cubit main window.
    \"\"\"
    app = QApplication.instance()
    if app is None:
        return None
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName() == 'claro':
            return w
    return None
```

## Usage pattern

```python
class CreateBrickDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # ... layout setup ...

# In the toolbar entry point:
if __name__ == '_coreform_cubit':
    claro = find_claro()  # may be None if run outside Cubit
    dlg = CreateBrickDialog(parent=claro)
    dlg.exec()
```

## Linter rule

`cubit-toolbar-missing-claro-parent` flags QDialog construction in
toolbar scripts that does not pass a `parent=` argument. This is a LOW
severity warning (not always needed for modal transient dialogs).

## Why objectName() and not windowTitle()

- `objectName()` is set programmatically by Cubit itself and is stable
  across versions and locales.
- `windowTitle()` can change (e.g., when a file is loaded Cubit may
  prepend the filename) and may be localized.

Always filter by `objectName() == 'claro'`, never by title.
"""

TOOLBAR_LLM_PROMPT_TEMPLATE = """
# LLM Prompt Template for Generating Cubit Toolbar Dialogs

Large language models (GPT-4, Claude, Gemini) generate PySide6 UI code
reliably when given the right preamble. The weak spot is the *Cubit*
parts - commands, APIs, idioms - because those are under-represented in
training data. Use this template and then have the LLM fill in the
dialog body.

## Canonical prompt template

```
Create a PySide6 dialog that runs inside Coreform Cubit.

HARD CONSTRAINTS (the Cubit embedded Python has quirks):
1. Begin the script with `#!python` as the first line.
2. DO NOT `import cubit` - the `cubit` namespace is already imported by
   the host.
3. When run inside Cubit, `__name__` is `'_coreform_cubit'` (not
   `'__main__'`). Guard the entry point with that check.
4. Cubit's Python importer DOES NOT support parenthesized multi-line
   imports. Use backslash continuations instead:
     from PySide6.QtWidgets import QDialog, QLabel, \\
         QLineEdit, QDialogButtonBox
5. Find the Cubit main window (QMainWindow objectName 'claro') and use
   it as the dialog parent to prevent the dialog from falling behind
   the main window.
6. Wrap float(line_edit.text()) conversions in try/except and show a
   QMessageBox on failure.
7. Drive Cubit with `cubit.cmd("...")` using f-strings (not an object
   API). This is the tested wire.

DIALOG SPEC:
- Title: "<title>"
- Row 1: QLabel "<label1>" + QLineEdit (default "<default1>")
- Row 2: QLabel "<label2>" + QLineEdit (default "<default2>")
- Buttons: OK / Cancel (QDialogButtonBox)
- On OK: read both line edits as floats, then execute:
    create brick x {v1}
    create cylinder radius {v2} height {v1 * 1.5}
    subtract volume 2 from volume 1

Please emit a complete script ready to paste into Cubit's custom
toolbar editor as a Python-script button.
```

## What the LLM tends to get wrong

| Mistake                              | Fix                                         |
|--------------------------------------|---------------------------------------------|
| `import cubit` at top                | Delete - already imported                    |
| `if __name__ == '__main__'`          | Replace with `'_coreform_cubit'`            |
| Parenthesized multi-line imports     | Use `\\` continuation                        |
| Dialog with no parent                | Add `find_claro()` + `parent=`              |
| `float()` without try/except         | Wrap in try/except + QMessageBox            |
| Calling `cubit.Volume(...)` objects  | Prefer `cubit.cmd("create volume ...")`     |

## Post-LLM checklist

1. Run `lint_cubit_script` on the output - catches most of the above.
2. Paste into the custom toolbar editor as an inline Python button.
3. Test with one parameter set that is known-good.
4. Export the toolbar + script as a `.tar.gz` once it works.
"""

TOOLBAR_PACKAGING = """
# Toolbar Packaging & Distribution

Once a toolbar works, package it for reuse and team sharing. Coreform
publishes two reference toolbars that show the exact layout.

## Reference toolbars

- `coreformLLC/cubit-toolbars` (GitHub): DAGMC + tire cross-section.
- Browse the layout before rolling your own.

## Directory layout

```
my_toolbar/
  toolbar.xml                  # toolbar definition (button names, icons,
                               # script paths, action types)
  scripts/
    create_brick.py            # Python script buttons
    tire_geometry.py
    cubit_util.py              # shared utility module (find_claro etc.)
  resources/
    icons/
      create_brick.png
      tire.png
  package-win.ps1              # Windows PowerShell packaging script
  package-linux.sh             # Linux bash packaging script
  README.md
```

## Packaging step

`package-win.ps1` / `package-linux.sh` simply `tar czf my_toolbar.tgz ...`
over the directory. No magic - the `.tar.gz` is a plain archive that
Cubit's toolbar importer unpacks under a namespace derived from the
archive filename.

## Import in Cubit

1. Open Cubit.
2. Right-click on the toolbar panel (where toolbars are listed).
3. Select "Import".
4. Browse to the `.tar.gz`.
5. Pick an extraction directory (e.g., a `plugins/` subfolder of your
   Cubit user home). Different toolbars land under their own namespace
   so multiple toolbars can coexist.
6. Enable the toolbar via Tools -> Custom Toolbar Editor.

## CI-friendly layout

If you want to build the toolbar archive in CI (so every push produces
a downloadable `.tar.gz`):

- Put `toolbar.xml` + `scripts/` + `resources/` at the root of a GitHub
  repo.
- Use a GitHub Actions workflow that runs `tar czf toolbar-{sha}.tgz
  toolbar.xml scripts resources` and attaches the artifact to the
  release.
- Users download the archive and import it the same way as the manual
  workflow.

## Sharing across a team

- Check the archive into a shared S: drive or Artifactory.
- Each user imports once; updates require re-importing (Cubit replaces
  the existing namespace).
- For lab-wide distribution on the 100goki farm, piggyback on the
  existing `deploy_radia_mcp_all_users.ps1` flow (see memory
  `project_mcp_multi_user_deploy.md`).
"""

TOOLBAR_ENCAPSULATION_PATTERN = """
# Encapsulation Pattern: Multi-Step Workflow as One Button

The biggest win of custom toolbars is turning a 10-step Cubit workflow
into a 3-click wizard. Coreform's tire cross-section toolbar is the
canonical example.

## Case study: tire cross-section toolbar

The user starts with just curves (no surface, no solid). The raw Cubit
workflow to mesh this is:

1. Import curves
2. Identify smallest curve (to choose merge tolerance)
3. Create a bounding surface
4. Imprint all curves onto the surface
5. Separate the resulting sub-surfaces
6. Delete the extraneous outer surface
7. Merge so interior edges are watertight
8. Identify which sub-surfaces are tread vs. base vs. ply (geometric
   heuristics)
9. Create blocks per material region
10. Set mesh size + mesh
11. Export

The toolbar compresses this to *three* buttons:

| Button | What it encapsulates                                            |
|--------|-----------------------------------------------------------------|
| 1. "Prepare geometry" | Steps 2-7 (auto-detects merge tolerance) |
| 2. "Assign blocks"    | Step 8-9 (auto-identifies tread/base)    |
| 3. "Mesh + export"    | Step 10-11 (default size, forced overwrite)|

## Design principles

1. **Auto-detect defaults** - scan geometry for smallest curve, suggest
   half that as merge tolerance. Offer a dialog to override, but ship
   a sensible default.
2. **One button = one logical unit of user intent**, not one Cubit
   command. "Prepare geometry" is intent; "imprint curves onto surface"
   is a Cubit command that serves that intent.
3. **Show the state after each click** - run `draw volume all with
   is_merged` or similar after a prepare step, so the user sees
   visual confirmation.
4. **Graceful degradation** - if auto-detection fails, pop the
   equivalent Cubit command-panel dialog instead of erroring out.
5. **Explicit buttons for the reversal path** - "Reset / undo last
   step" buttons matter when users are iterating on parameters.

## Skeleton script

```python
#!python
# Encapsulated step: prepare geometry for tire cross-section

from PySide6.QtWidgets import QDialog, QLabel, QLineEdit, QVBoxLayout, \\
    QDialogButtonBox, QMessageBox, QApplication, QMainWindow

def find_claro():
    app = QApplication.instance()
    if app is None:
        return None
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName() == 'claro':
            return w
    return None

def smallest_curve_length():
    \"\"\"Query Cubit for the shortest curve in the current session.\"\"\"
    ids = cubit.parse_cubit_list("curve", "all")
    lengths = [cubit.curve(cid).length() for cid in ids]
    if not lengths:
        return None
    return min(lengths), ids[lengths.index(min(lengths))]

class PrepareDialog(QDialog):
    def __init__(self, default_tol, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Prepare Geometry")
        self.tol_edit = QLineEdit(str(default_tol))
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f"Merge tolerance (default {default_tol}):"))
        lay.addWidget(self.tol_edit)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

def prepare_geometry():
    result = smallest_curve_length()
    if result is None:
        QMessageBox.warning(None, "No geometry",
                            "Import curves before pressing Prepare.")
        return
    length, cid = result
    default_tol = length / 2.0

    claro = find_claro()
    dlg = PrepareDialog(default_tol, parent=claro)
    if dlg.exec() != QDialog.Accepted:
        return

    try:
        tol = float(dlg.tol_edit.text())
    except ValueError as exc:
        QMessageBox.warning(claro, "Invalid input", str(exc))
        return

    cubit.cmd(f"merge tolerance {tol}")
    cubit.cmd("create surface curve all")
    cubit.cmd("imprint all")
    cubit.cmd("merge all")
    cubit.cmd("delete surface with area > 1e6")  # drop outer bounding surf

if __name__ == '_coreform_cubit':
    prepare_geometry()
```
"""

TOOLBAR_QUALITY_CONVERGENCE = """
# Quality Convergence Loop Pattern

A recurring question from Coreform webinars: "can I auto-mesh, query
quality, and refine until converged?" Yes - this is a natural fit for a
toolbar button. Here is the canonical shape.

## Loop structure

```python
#!python
# Auto-refine mesh until all elements pass a shape threshold or a
# refinement budget is exhausted.

from PySide6.QtWidgets import QMessageBox

SHAPE_THRESHOLD = 0.2      # Cubit's `shape` metric, 0.2 is the
                           # conventional "usable" cutoff
MAX_ITERATIONS = 5         # hard budget - never run forever
REFINE_FACTOR = 0.7        # size shrink per iteration
INITIAL_SIZE = 1.0

def min_shape_quality():
    \"\"\"Return the minimum `shape` quality across all tet elements.\"\"\"
    # Evaluate quality on all volumes; Cubit caches result
    cubit.cmd("quality volume all shape global")
    # Use the Python API to read element-wise quality
    tets = cubit.parse_cubit_list("tet", "all")
    if not tets:
        return None
    qs = [cubit.get_quality_value("tet", t, "shape") for t in tets]
    return min(qs), qs.index(min(qs))

def auto_refine(volume_id):
    size = INITIAL_SIZE
    cubit.cmd(f"volume {volume_id} scheme tetmesh")
    for it in range(MAX_ITERATIONS):
        cubit.cmd(f"volume {volume_id} size {size}")
        cubit.cmd(f"mesh volume {volume_id}")
        result = min_shape_quality()
        if result is None:
            return False, "No elements"
        q, bad_tet_idx = result
        if q >= SHAPE_THRESHOLD:
            return True, f"Converged at iter {it+1}, min shape = {q:.3f}"
        # Shrink size and remesh
        size *= REFINE_FACTOR
        cubit.cmd(f"delete mesh volume {volume_id} propagate")
    return False, (
        f"Budget exhausted ({MAX_ITERATIONS} iters), "
        f"last min shape = {q:.3f} < threshold {SHAPE_THRESHOLD}"
    )

if __name__ == '_coreform_cubit':
    ok, msg = auto_refine(volume_id=1)
    QMessageBox.information(None, "Auto-refine", msg)
```

## Design rules

1. **Always budget the loop.** A hard iteration cap (MAX_ITERATIONS)
   prevents runaway refinement on pathological geometry.
2. **Check for trivial failure first.** `min_shape_quality` should
   return None if there are zero elements rather than crashing.
3. **Scope the refinement.** Shrinking `size` on the whole volume is
   coarse; for production, refine only in the neighborhood of the bad
   element: `refine tet {bad_tet_idx} numsplit 1`.
4. **Report both success and failure.** Use QMessageBox so the user
   sees the final metric, not just "done."
5. **Delete the mesh between iterations.** `delete mesh volume N
   propagate` is the right reset; plain `delete mesh` does not clear
   attached settings.

## Quality metrics quick reference

| Metric           | Good (tet) | Good (hex) | Notes                       |
|------------------|-----------|-----------|-----------------------------|
| `shape`          | >= 0.2    | >= 0.3    | Cubit default, meta-quality |
| `scaled jacobian`| >= 0.2    | >= 0.5    | Rigor for Jacobian ratio    |
| `aspect ratio`   | <= 3      | <= 3      | Lower is better             |
| `skew`           | <= 0.5    | <= 0.5    | 0 is perfect                |
| `distortion`     | >= 0.6    | >= 0.6    | Higher is better            |

`shape` is the right first choice because it rolls several metrics into
one heuristic. Escalate to `scaled_jacobian` if the solver (Abaqus,
Ansys, radia/NGSolve) reports Jacobian warnings.
"""

TOOLBAR_TROUBLESHOOTING = """
# Troubleshooting In-Cubit Toolbar Scripts

## Symptom: "The toolbar button does nothing"

- Is `__name__ == '_coreform_cubit'` wrapping the entry call? A script
  whose whole body is under `if __name__ == '__main__':` never runs.
- Did Cubit open a dialog that fell behind the main window? Fix with
  `find_claro()` + `parent=claro`.

## Symptom: "ImportError on PySide6 imports"

- Parenthesized multi-line imports break Cubit's importer. Switch to
  backslash continuations or single-line imports.
- Radia targets Coreform Cubit 2025.12+. Do not add a PyQt5 fallback or
  patch older Cubit bundles; upgrade Cubit and re-run `cubit-plugin-install`.

## Symptom: "cubit.cmd raises 'not initialized'"

- If running outside Cubit for unit tests, you need
  `cubit.init(['cubit', '-nojournal', '-batch'])`. Do NOT include that
  line inside the toolbar script - Cubit is already initialized when
  the script runs there. Gate it on
  `if __name__ != '_coreform_cubit'`.

## Symptom: "Float conversion crashes the Cubit terminal"

- `float(line_edit.text())` without try/except will spew a traceback to
  the terminal and not tell the user. Always wrap in try/except and
  show a QMessageBox.

## Symptom: "Toolbar import fails to load scripts"

- Check the `.tar.gz` layout - scripts must be under `scripts/`, icons
  under `resources/`, toolbar config as `toolbar.xml` at the archive
  root. Extra top-level directories from `tar` (e.g., tarring the
  parent folder by mistake) break the importer.

## Symptom: "LLM-generated dialog works once but breaks after code
   regeneration"

- LLMs often "fix" a working script by adding `import cubit` or
  changing `__name__` checks. Lock the working version in your
  toolbar's `scripts/` folder and only regenerate the dialog portion,
  not the Cubit glue.

## Symptom: "Dialog becomes unresponsive during long Cubit command"

- Cubit commands block the Qt event loop. For long operations, either
  (a) pre-estimate duration and show a `QProgressDialog`, or (b) break
  the work into smaller `cubit.cmd` calls with `QApplication.
  processEvents()` between them. Do NOT spawn a QThread - Cubit's
  Python binding is not thread-safe.
"""

# ---- Topic registry ----------------------------------------------------

def get_toolbar_documentation(topic: str = "overview") -> str:
    """Return in-Cubit custom toolbar documentation by topic."""
    topics = {
        "overview": TOOLBAR_OVERVIEW,
        "python_conventions": TOOLBAR_PYTHON_SCRIPT_CONVENTIONS,
        "conventions": TOOLBAR_PYTHON_SCRIPT_CONVENTIONS,  # alias
        "shebang": TOOLBAR_PYTHON_SCRIPT_CONVENTIONS,       # alias
        "import_gotcha": TOOLBAR_PYTHON_SCRIPT_CONVENTIONS, # alias
        "claro": TOOLBAR_CLARO_PARENT_HELPER,
        "parent": TOOLBAR_CLARO_PARENT_HELPER,              # alias
        "find_claro": TOOLBAR_CLARO_PARENT_HELPER,          # alias
        "llm_prompt": TOOLBAR_LLM_PROMPT_TEMPLATE,
        "prompt_template": TOOLBAR_LLM_PROMPT_TEMPLATE,     # alias
        "dialog_prompt": TOOLBAR_LLM_PROMPT_TEMPLATE,       # alias
        "packaging": TOOLBAR_PACKAGING,
        "distribution": TOOLBAR_PACKAGING,                  # alias
        "tar_gz": TOOLBAR_PACKAGING,                        # alias
        "encapsulation": TOOLBAR_ENCAPSULATION_PATTERN,
        "workflow_pattern": TOOLBAR_ENCAPSULATION_PATTERN,  # alias
        "tire_example": TOOLBAR_ENCAPSULATION_PATTERN,      # alias
        "quality_convergence": TOOLBAR_QUALITY_CONVERGENCE,
        "auto_refine": TOOLBAR_QUALITY_CONVERGENCE,         # alias
        "convergence_loop": TOOLBAR_QUALITY_CONVERGENCE,    # alias
        "troubleshooting": TOOLBAR_TROUBLESHOOTING,
        "debug": TOOLBAR_TROUBLESHOOTING,                   # alias
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join([
            TOOLBAR_OVERVIEW,
            TOOLBAR_PYTHON_SCRIPT_CONVENTIONS,
            TOOLBAR_CLARO_PARENT_HELPER,
            TOOLBAR_LLM_PROMPT_TEMPLATE,
            TOOLBAR_PACKAGING,
            TOOLBAR_ENCAPSULATION_PATTERN,
            TOOLBAR_QUALITY_CONVERGENCE,
            TOOLBAR_TROUBLESHOOTING,
        ])
    if topic in topics:
        return topics[topic]
    return (
        f"Unknown topic: '{topic}'. "
        f"Available: all, overview, python_conventions, claro, "
        f"llm_prompt, packaging, encapsulation, quality_convergence, "
        f"troubleshooting. (Aliases accepted: import_gotcha, "
        f"find_claro, prompt_template, tar_gz, tire_example, "
        f"auto_refine, debug.)"
    )


# ---- Skeleton generators ------------------------------------------------

def generate_toolbar_skeleton(name: str, buttons: list[dict]) -> dict:
    """Generate a complete custom-toolbar skeleton (dict of relative path
    -> file content).

    Args:
        name: toolbar name, e.g. "radia_ih_helpers"
        buttons: list of dicts with keys:
            - 'label': button label shown in Cubit toolbar
            - 'kind': 'python' | 'cubit_script'
            - 'script_basename': stem of the script filename (no ext)
            - 'body': initial script body (without shebang/boilerplate)

    Returns:
        dict mapping relative paths -> file content (str). Caller writes
        these to disk and then runs the platform packaging script.
    """
    import textwrap

    safe_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)
    files: dict[str, str] = {}

    # toolbar.xml (minimal, one toolbar with the given buttons)
    button_entries = []
    for b in buttons:
        kind = b.get("kind", "python")
        label = b.get("label", "Untitled")
        basename = b.get("script_basename", "button")
        if kind == "python":
            script_rel = f"scripts/{basename}.py"
            button_entries.append(
                f'    <button label="{label}" type="python" '
                f'script="{script_rel}" icon="resources/icons/{basename}.png"/>'
            )
        else:
            script_rel = f"scripts/{basename}.jou"
            button_entries.append(
                f'    <button label="{label}" type="cubit_script" '
                f'script="{script_rel}" icon="resources/icons/{basename}.png"/>'
            )

    files["toolbar.xml"] = (
        f'<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<toolbar name="{safe_name}">\n'
        + "\n".join(button_entries)
        + "\n</toolbar>\n"
    )

    # Shared utility with find_claro()
    files["scripts/cubit_util.py"] = textwrap.dedent('''\
        #!python
        """Shared utilities for the toolbar. Import with:

            from cubit_util import find_claro
        """

        from PySide6.QtWidgets import QApplication, QMainWindow


        def find_claro():
            """Return Cubit's main window (QMainWindow named 'claro')."""
            app = QApplication.instance()
            if app is None:
                return None
            for w in app.topLevelWidgets():
                if isinstance(w, QMainWindow) and w.objectName() == 'claro':
                    return w
            return None
    ''')

    # Each Python button: shebang + claro + try/except wrapper
    for b in buttons:
        kind = b.get("kind", "python")
        basename = b.get("script_basename", "button")
        label = b.get("label", "Untitled")
        body = b.get("body", f"cubit.cmd('echo \"hello from {label}\"')")
        if kind == "python":
            indent_body = "\n".join(
                " " * 12 + line for line in body.splitlines()
            )
            header = (
                "#!python\n"
                f'"""Auto-generated skeleton for toolbar button: {label}\n'
                "\n"
                "- __name__ == '_coreform_cubit' when run from the toolbar.\n"
                "- `cubit` is already imported by the host; do not re-import.\n"
                "- Cubit's in-process importer dislikes parenthesized imports;\n"
                "  use backslash continuations.\n"
                '"""\n\n'
                "from PySide6.QtWidgets import QMessageBox, QDialog, QLabel, \\\n"
                "    QLineEdit, QVBoxLayout, QDialogButtonBox\n"
                "\n"
                "from cubit_util import find_claro\n"
                "\n\n"
                "def _run():\n"
                "    parent = find_claro()\n"
                "    try:\n"
            )
            footer = (
                "    except Exception as exc:\n"
                f'        QMessageBox.warning(parent, "{label} - error",\n'
                '                            f"{type(exc).__name__}: {exc}")\n'
                "\n\n"
                "if __name__ == '_coreform_cubit':\n"
                "    _run()\n"
            )
            files[f"scripts/{basename}.py"] = header + indent_body + "\n" + footer
        else:  # cubit_script
            files[f"scripts/{basename}.jou"] = (
                f"## Auto-generated journal for toolbar button: {label}\n"
                + body
                + "\n"
            )

    # README
    files["README.md"] = textwrap.dedent(f'''\
        # {safe_name}

        Custom Cubit toolbar scaffold generated by `cubit_scaffold_toolbar`.

        ## Layout

        ```
        {safe_name}/
          toolbar.xml               # Toolbar definition
          scripts/                  # Python + journal scripts
            cubit_util.py           # find_claro() helper
          resources/icons/          # Button icons (add PNGs here)
          README.md
        ```

        ## Package and install

        ```bash
        # from the parent directory of this folder:
        tar czf {safe_name}.tgz {safe_name}/
        # Open Cubit -> right-click toolbar panel -> Import -> {safe_name}.tgz
        ```

        ## Buttons

        ''' + "\n".join(f"- **{b['label']}** (`{b.get('kind','python')}`): "
                        f"`scripts/{b.get('script_basename','button')}"
                        f"{'.py' if b.get('kind','python')=='python' else '.jou'}`"
                        for b in buttons) + "\n")

    # Windows packaging script
    files["package-win.ps1"] = textwrap.dedent(f'''\
        # PowerShell 7+ - produces {safe_name}.tgz in the parent directory.
        $ErrorActionPreference = 'Stop'
        $here = Split-Path -Parent $MyInvocation.MyCommand.Path
        $parent = Split-Path -Parent $here
        $name = Split-Path -Leaf $here
        Push-Location $parent
        try {{
            tar czf "$name.tgz" $name
            Write-Host "Wrote $parent\\$name.tgz"
        }} finally {{
            Pop-Location
        }}
    ''')

    # Linux packaging script
    files["package-linux.sh"] = textwrap.dedent(f'''\
        #!/usr/bin/env bash
        set -euo pipefail
        here="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
        parent="$(dirname "$here")"
        name="$(basename "$here")"
        cd "$parent"
        tar czf "$name.tgz" "$name"
        echo "Wrote $parent/$name.tgz"
    ''')

    return files


def generate_dialog_skeleton(title: str, fields: list[dict],
                             on_accept_commands: list[str]) -> str:
    """Generate a complete PySide6 dialog script for a Cubit toolbar
    button.

    Args:
        title: dialog window title.
        fields: list of dicts, each with keys:
            - 'name': attribute name for the QLineEdit (e.g., 'size').
            - 'label': text label shown next to the field.
            - 'default': default text in the QLineEdit.
            - 'kind': 'float' | 'int' | 'text' (controls validation).
        on_accept_commands: list of `cubit.cmd` argument strings. Each
            may use `{<field_name>}` placeholders that will be f-string
            substituted at runtime with the validated field values.

    Returns:
        Complete Python script as a string, including shebang, claro
        helper, dialog class, and entry point.
    """
    # Dialog class name (camel-case)
    dlg_cls = "".join(w.capitalize() for w in title.split()) + "Dialog"

    # Field declarations inside __init__ (8-space indent)
    init_field_lines = []
    for f in fields:
        name = f["name"]
        label = f.get("label", name)
        default = f.get("default", "")
        init_field_lines.append(f'        self.{name}_edit = QLineEdit("{default}")')
        init_field_lines.append(f'        lay.addWidget(QLabel("{label}:"))')
        init_field_lines.append(f'        lay.addWidget(self.{name}_edit)')

    # Field reads inside _run try block (8-space indent: the 'try:' is
    # at 4 spaces, contents at 8)
    read_field_lines = []
    for f in fields:
        name = f["name"]
        kind = f.get("kind", "text")
        if kind == "float":
            read_field_lines.append(f'        {name} = float(dlg.{name}_edit.text())')
        elif kind == "int":
            read_field_lines.append(f'        {name} = int(dlg.{name}_edit.text())')
        else:
            read_field_lines.append(f'        {name} = dlg.{name}_edit.text()')

    # Cubit command calls inside _run, after try/except (4-space indent)
    cmd_lines = []
    for cmd in on_accept_commands:
        if "{" in cmd and "}" in cmd:
            cmd_lines.append(f'    cubit.cmd(f"{cmd}")')
        else:
            cmd_lines.append(f'    cubit.cmd("{cmd}")')

    parts = [
        "#!python",
        f'"""Auto-generated Cubit toolbar dialog: {title}',
        "",
        "Generated by radia_mcp.cubit.generate_dialog_skeleton.",
        '"""',
        "",
        "from PySide6.QtWidgets import QApplication, QDialog, QDialogButtonBox, \\",
        "    QLabel, QLineEdit, QMainWindow, QMessageBox, QVBoxLayout",
        "",
        "",
        "def find_claro():",
        "    app = QApplication.instance()",
        "    if app is None:",
        "        return None",
        "    for w in app.topLevelWidgets():",
        "        if isinstance(w, QMainWindow) and w.objectName() == 'claro':",
        "            return w",
        "    return None",
        "",
        "",
        f"class {dlg_cls}(QDialog):",
        "    def __init__(self, parent=None):",
        "        super().__init__(parent)",
        f'        self.setWindowTitle("{title}")',
        "        lay = QVBoxLayout(self)",
    ]
    parts.extend(init_field_lines)
    parts.extend([
        "        bb = QDialogButtonBox(",
        "            QDialogButtonBox.Ok | QDialogButtonBox.Cancel",
        "        )",
        "        bb.accepted.connect(self.accept)",
        "        bb.rejected.connect(self.reject)",
        "        lay.addWidget(bb)",
        "",
        "",
        "def _run():",
        "    parent = find_claro()",
        f"    dlg = {dlg_cls}(parent=parent)",
        "    if dlg.exec() != QDialog.Accepted:",
        "        return",
        "    try:",
    ])
    parts.extend(read_field_lines)
    parts.extend([
        "    except ValueError as exc:",
        '        QMessageBox.warning(parent, "Invalid input",',
        '                            f"{type(exc).__name__}: {exc}")',
        "        return",
    ])
    parts.extend(cmd_lines)
    parts.extend([
        "",
        "",
        "if __name__ == '_coreform_cubit':",
        "    _run()",
        "",
    ])
    return "\n".join(parts)
