"""
Cubit-bypass standalone launch of the four Radia-NGSolve panels —
how to run radia-ih / radia-em / radia-pcb / radia-heat as plain
PySide6 desktop apps without Cubit installed.

Read this when:
* A user has a `.vol` mesh from any source (Cubit, Netgen-OCC,
  build123d, etc.) and wants to run an analysis panel without paying
  the Coreform Cubit licence.
* Setting up a non-Cubit lab seat (LAB-side dev clone, mdx, a
  collaborator's machine).
* Diagnosing "the panel I see in Cubit's Solve menu — can I run it
  outside Cubit?"  (Yes, via the standalone .exe entry-points.)

The MCP server exposes this via standalone_panels(topic=...). Topics:
quick_start, four_panels, vol_sources, vs_cubit, ih_methods,
troubleshooting.
"""

STANDALONE_PANELS = """\
# Radia-NGSolve standalone panel launch (Cubit-bypass)

The four PySide6 panels (radia-ih, radia-em, radia-pcb, radia-heat)
are designed to run **as plain desktop apps without Cubit**.  Each
ships as a console entry-point installed by
`pip install radia[cubit,gui]` (or `radia[gui]` if you don't need the
Cubit plugin extras).  The panels read a `.vol` mesh that you supply
— Cubit is only one of many `.vol` sources.

Topics: quick_start, four_panels, vol_sources, vs_cubit, ih_methods,
        troubleshooting

============================================================
## quick_start — Kubota recipe: 100% Cubit-free panel launch
============================================================

End-to-end recipe for a lab member (e.g. Kubota) who wants to run an
analysis panel **without Cubit installed at all** — no Coreform Cubit
licence, no `cubit-plugin-install`, no Cubit process running in the
background.  Cubit literally does not appear anywhere in this flow.

```cmd
:: 1. install — `[gui]` extras only.  Note: NOT `[cubit,gui]`.
::    `[cubit]` would pull cubit-mesh-export which only matters if
::    you ALSO want the in-Cubit Solve menu launcher.  For a
::    Cubit-free workflow, `[gui]` alone is enough.
pip install --upgrade "radia[gui]" radia-mcp

:: 2. confirm the four .exe launchers are on PATH.
where radia-ih radia-em radia-pcb radia-heat

:: 3. launch the IH panel pre-pointed at your .vol.
::    The panel pops up directly as a desktop app.
radia-ih C:\\path\\to\\model.vol

:: (optional) launch with no argument and use the in-panel
::            Working folder + per-method Browse fields instead.
radia-ih
```

Where does `model.vol` come from?  See the `vol_sources` topic — the
short answer is **NGSolve OCC standalone** (Python-only, no Cubit
licence; takes a STEP file in, writes a `.vol` out).  Once the `.vol`
exists, steps 2 and 3 above are the entire panel workflow.

Inside the IH panel:

| # | Action | Notes |
|---|--------|-------|
| 1 | **Method** dropdown — pick one of 9 | Six EM methods plus three Thermal methods integrated into `radia-ih`. |
| 2 | Browse the visible geometry rows | PEEC modes use coil STEP. BEM-A modes use a pre-meshed coil `.vol` with `source` / `sink`. Workpiece modes use the panel-owned Workpiece `.vol` row. Thermal modes use workpiece `.vol` + qsurf `.sol` + EM `.vol`. |
| 3 | **Frequency [Hz]** + **Coil current [A]** | Operating point for skin depth / energy / dissipation. |
| 4 | **Coil material** + **Workpiece material** | Pick from preset (Cu / Steel mu_r=100 / etc.) or "Custom" + type sigma / mu_r. |
| 5 | **Solver** | Dense LU (small) / HACApK (large) — affects assembly + memory. |
| 6 | Click **Run** | A Layer 4 subprocess (`calc_inductance.py` / `calc_fem_kelvin.py` etc.) is launched.  Output appears in the lower text area as JSON. |
| 7 | (optional) **Open GMSH** | For modes that emit `.msh`, opens the workpiece + coil in a GMSH viewer. |

The Cubit launcher is unnecessary for steps 1-7 above.  Cubit is
only required when you want **the panel to drive a `.jou` file +
mesh export from a `.cub5`** — that workflow is exposed via Cubit's
`Solve -> Radia-NGSolve` menu, NOT via the standalone .exe.  Once a
`.vol` exists, the standalone launcher and the Cubit launcher are
identical (both call the same Layer 4 `calc_*.py` subprocess).

The top row is a **Working folder**.  Browse-row paths are displayed
relative to it when possible, then resolved back to absolute paths
before subprocess launch.

============================================================
## four_panels — what each launcher does
============================================================

| Launcher       | Window class | Purpose                                                        | Required input             |
|----------------|--------------|----------------------------------------------------------------|----------------------------|
| `radia-em`     | `EMWindow`   | Electromagnet: Omega-reduced / A-Phi / MSC / Kelvin Benchmark | `.vol`                     |
| `radia-ih`     | `IHWindow`   | Induction Heating + Thermal: PEEC / BEM-A / FEM A-V + q_surf->T post | STEP (PEEC) + `.vol` (BEM/FEM/thermal) |
| `radia-pcb`    | `PCBWindow`  | PCB / FastHenry impedance sweep                                | FastHenry `.inp`           |
| `radia-motor`  | `MotorWindow`| Motor: transient (Lange-Henrotte-Hameyer) + lamination (Hollaus EM) | `.vol`                  |

All four are top-level PySide6 `QMainWindow` apps.  Launching them
spawns a single Python process (the Layer 3 panel); when the user
clicks **Run** that process spawns a Layer 4 `calc_*.py`
subprocess.  Cubit is not in either chain.

============================================================
## vol_sources — generate .vol without Cubit
============================================================

The standalone panels accept any Netgen `.vol` (text format) that
declares the labels each method expects.  For a **fully Cubit-free**
pipeline (no licence, no plugin, no Cubit process), use NGSolve OCC
or build123d directly.

### NGSolve OCC standalone (recommended Cubit-free path)

Pure Python, ships with `pip install ngsolve` (already a transitive
dependency of `radia[gui]`).  Reads a STEP, attaches face names for
the panel-required labels, writes a `.vol`:

```python
# make_vol_ih.py  --  example: IH workpiece with sibc sideset.
from netgen.occ import OCCGeometry
from netgen.meshing import MeshingParameters

geo = OCCGeometry("workpiece.step")
shape = geo.shape

# Tag faces.  Panel methods read these as boundary labels.
# For PEEC + BEM weak (radia-ih) the workpiece needs `sibc` on its
# whole outer surface.
for f in shape.faces:
    f.name = "sibc"

# Volumetric label (becomes the .vol material).
shape.solids.first.name = "workpiece"

ngmesh = OCCGeometry(shape).GenerateMesh(
    mp=MeshingParameters(maxh=0.005, curvaturesafety=1.0))
ngmesh.Save("workpiece.vol")
print("wrote workpiece.vol")
```

For FEM A-V full mode, the coil also needs `coil` / `source` / `sink`
labels — see `examples/induction_heating/` for templates.

### build123d → STEP → NGSolve OCC

For programmatic geometry (parametric coils, parametric workpieces),
build geometry in build123d, write STEP, mesh via NGSolve OCC.  See
`tests/bem/fixtures/make_gapped_torus_2port.py` for a working
example used by the BEM-A tests.

### From a colleague over W:\\ / e-mail

If a `.vol` already exists (NAS share, attachment, etc.) just pass
its path to the launcher.  Cubit's authorship of the `.vol` is
irrelevant — once written, a `.vol` is a plain text file the panel
reads directly.

### `.vol` label requirements per panel method

The `.vol` is text and can be inspected with any editor; the
required label sets are documented per-method in the panel's
`check_method_requirements()`:

| Panel method                          | Required volume materials | Required boundaries          |
|---------------------------------------|---------------------------|-------------------------------|
| PEEC inductance (vacuum)              | none                      | none                          |
| BEM-A inductance (vacuum)             | none                      | none                          |
| PEEC + BEM weak coupling              | none                      | `sibc`                        |
| BEM-A + BEM weak coupling             | none                      | `sibc`                        |
| PEEC coil + FEM Kelvin                | `kelvin`                  | `sibc`                        |
| FEM A-V full                          | `coil`, `kelvin`          | `source`, `sink`, `sibc`      |

Tip: launch the panel first with no `.vol`, pick the method, then
read the status line under the Method dropdown — it lists exactly
which labels are missing or present in your `.vol`.

============================================================
## vs_cubit — when (if ever) you need Cubit
============================================================

The standalone .exe launchers (`radia-ih` / `radia-em` / `radia-pcb`
/ `radia-heat`) are the **default and only path** that this knowledge
module documents.  Cubit is mentioned here only to clarify what is
NOT needed:

* **Cubit licence** — not needed for the standalone launchers.
* **`cubit-plugin-install`** — not needed.  Skip it entirely on a
  Cubit-free machine.
* **A running Cubit process** — never started.  The standalone
  panels are top-level PySide6 desktop apps.
* **`[cubit]` extras in pip install** — not needed; use `[gui]`
  alone.  `[cubit,gui]` only adds value for users who also want the
  in-Cubit Solve-menu launcher.

A separate launch path exists inside Coreform Cubit (`Solve ->
Radia-NGSolve`) for users who already use Cubit for mesh generation
and prefer to stay inside Cubit's GUI.  That path is documented in
`cubit_plugin_layers` (different topic, different audience).  For
the Kubota-style "open my .vol, run the analysis, close the window"
workflow this knowledge module is concerned with, no Cubit-side
component is involved.

Practical implication: the standalone-panel install on a developer
laptop / non-Cubit lab seat is a single `pip install` line
(`radia[gui]==X.Y.Z`) plus whatever generates the `.vol` (NGSolve
OCC standalone — see `vol_sources`).  Total install state: a
Python interpreter + the radia wheel + PySide6.

============================================================
## ih_methods — 9 IH methods on radia-ih (EM + integrated thermal)
============================================================

Since 4.63.0 the IH panel offers 9 methods.  The first six are EM
methods; the last three run the integrated thermal follow-up.  Two EM
pairs differ only in the coil solver (PEEC perimeter filaments vs
BEM-A surface RWG):

| Method                                | Coil          | Workpiece          | Inputs                  | Cost     |
|---------------------------------------|---------------|--------------------|-------------------------|----------|
| PEEC inductance (vacuum)              | PEEC filament | none               | STEP only               | ~30 s    |
| BEM-A inductance (vacuum)             | BEM-A surface | none               | coil `.vol`             | ~25 s    |
| PEEC + BEM weak coupling              | PEEC filament | scalar BEM-SIBC    | STEP + `.vol`           | ~3 min   |
| BEM-A + BEM weak coupling             | BEM-A surface | scalar BEM-SIBC    | coil `.vol` + wp `.vol` | ~3-5 min |
| PEEC coil + FEM wp (SIBC) + Kelvin    | PEEC filament | volumetric FEM SIBC | STEP + `.vol` (Kelvin) | ~4-8 min |
| Full simulation (FEM A-V + wp SIBC)   | volumetric    | volumetric FEM SIBC | `.vol` (full labels)    | ~1-7 min |
| Thermal: 3D static                    | none          | heat equation       | wp `.vol` + qsurf `.sol` + EM `.vol` | varies |
| Thermal: 3D + rotation                | none          | rotating heat       | wp `.vol` + qsurf `.sol` + EM `.vol` | varies |
| Thermal: 2D axisymmetric              | none          | axisym heat         | axisym wp `.vol` + qsurf `.sol` + EM `.vol` | fast |

PEEC vacuum mode needs only a coil STEP; BEM-A vacuum mode needs a
pre-meshed coil `.vol` with `source` and `sink` terminal labels.  Both
skip the workpiece `.vol`.

For weak-coupled and FEM modes, the `.vol` is required and must
have the right labels (`sibc` boundary for BEM-SIBC; `coil` /
`source` / `sink` / `sibc` / `kelvin` for FEM A-V full).

File pickers resolve against the top-level **Working folder**.  The
workpiece mesh is the panel-owned `wp_vol` row, not the hidden legacy
`AnalysisWindow._vol_edit` shim.

============================================================
## troubleshooting — common errors when launching standalone
============================================================

### `ModuleNotFoundError: No module named 'PySide6'`

You installed `radia` or `radia[cubit]` without the `[gui]` extra.
PySide6 is in `[gui]`.  Fix:

```cmd
pip install --upgrade "radia[cubit,gui]==X.Y.Z"
```

### `radia-ih: command not found` (or "is not recognized" on Windows)

The console entry-point is in `<Python>/Scripts/` but that directory
is not on PATH.

* Windows: `C:\\Program Files\\Python312\\Scripts\\` should be on
  PATH; if not, run the launcher with its full path.
* Use `where radia-ih` (Windows) / `which radia-ih` (Linux/macOS)
  to confirm.

### Panel opens but `.vol` Browse field is greyed out

You picked a method that doesn't use a `.vol` (e.g. PEEC inductance
in IH, or PCB).  This is correct behaviour — those modes don't need
a workpiece mesh.

### "Boundary 'sibc' not found in .vol" / "Material 'coil' not found"

The `.vol` you supplied is missing labels the chosen method
requires.  The status line under the Method dropdown lists exactly
which labels are missing.  Either re-export the `.vol` with the
right block / sideset names (Cubit) or pick a method whose
requirements your `.vol` already satisfies.

### Run button produces no output and no error

The Layer 4 subprocess crashed before emitting JSON.  Check the
panel's lower text area (called "Output") for stderr.  Most common
cause: a NumPy / NGSolve import error from a stale install — run
`python -c "import radia, ngsolve; print(radia.__version__,
ngsolve.__version__)"` from the same Python the panel is using.

### `radia-heat` is missing from `where` / `which` output

`radia-heat` was added in 4.25.1.  `pip install --upgrade
radia[gui]==4.25.1` (or later) registers it.  4.25.0 had only three
launchers (radia-em, radia-ih, radia-pcb).
"""


_TOPICS = (
    "quick_start", "four_panels", "vol_sources", "vs_cubit",
    "ih_methods", "troubleshooting",
)


def get_standalone_panels_documentation(topic: str = "") -> str:
    """Return the standalone-panel knowledge.

    Args:
        topic: Empty for the full document, or one of the entries in
               ``_TOPICS`` above for a single section.

    Returns:
        Markdown string.  When ``topic`` is empty the entire document
        is returned (~10 kB); when a topic is named only that section
        is sliced out (between its `## topic ` header and the next
        topic header).  The leading `=========` separator is included
        on the slice for visual continuity.
    """
    if not topic:
        return STANDALONE_PANELS

    topic = topic.strip()
    if topic not in _TOPICS:
        return (f"Unknown topic {topic!r}.  Valid topics: "
                + ", ".join(_TOPICS))

    # Locate every "\n## <one-of-topic> " header offset so we can
    # slice between consecutive topics.
    headers = []
    pos = 0
    while True:
        nxt = STANDALONE_PANELS.find("\n## ", pos)
        if nxt < 0:
            break
        line_end = STANDALONE_PANELS.find("\n", nxt + 1)
        line = STANDALONE_PANELS[nxt + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, nxt + 1))
                break
        pos = nxt + 1

    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."

    # Include the `=========` separator immediately above the topic
    # header so the rendered slice keeps its visual frame.
    section_start = req_starts[0]
    sep_above = STANDALONE_PANELS.rfind(
        "============================================================",
        0, section_start)
    if sep_above >= 0 and sep_above > section_start - 80:
        section_start = sep_above

    # End at the next topic's separator-above (or document end).
    last_req = req_starts[-1]
    section_end = len(STANDALONE_PANELS)
    for kw, off in headers:
        if kw != topic and off > last_req:
            sep = STANDALONE_PANELS.rfind(
                "============================================================",
                0, off)
            section_end = sep if sep > 0 else off
            break

    return STANDALONE_PANELS[section_start:section_end].rstrip() + "\n"
