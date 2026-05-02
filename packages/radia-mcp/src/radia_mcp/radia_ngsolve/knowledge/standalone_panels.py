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
## quick_start — Kubota recipe: run radia-ih on an existing .vol (no Cubit)
============================================================

End-to-end recipe for a lab member (e.g. Kubota) who has a `.vol`
mesh and wants to compute IH coil + workpiece quantities without
opening Cubit:

```cmd
:: 1. install (one-time per machine, machine-wide on 100号機).
::    [cubit,gui] extras pull cubit-mesh-export AND PySide6.
::    [gui] alone is enough if you don't need the Cubit plugin.
pip install --upgrade "radia[cubit,gui]==4.25.1" radia-mcp==0.36.7

:: 2. confirm the four .exe launchers are on PATH.
where radia-ih radia-em radia-pcb radia-heat

:: 3. launch the IH panel pre-pointed at your .vol.
::    No Cubit process needed; the panel pops up directly.
radia-ih C:\\path\\to\\model.vol

:: (optional) launch with no argument and use the in-panel
::            "Model (.vol):  [Browse...]" field instead.
radia-ih
```

Inside the IH panel:

| # | Action | Notes |
|---|--------|-------|
| 1 | **Method** dropdown — pick one of 6 | PEEC inductance / BEM-A inductance / PEEC + BEM weak / BEM-A + BEM weak / PEEC + FEM Kelvin / FEM A-V full |
| 2 | **STEP or JOU:**  [Browse...] (if PEEC variants) | The coil STEP file.  PEEC + BEM and BEM-A + BEM weak ALSO need a `.vol` for the workpiece — the launcher arg above. |
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

============================================================
## four_panels — what each launcher does
============================================================

| Launcher       | Window class | Purpose                                                        | Required input             |
|----------------|--------------|----------------------------------------------------------------|----------------------------|
| `radia-em`     | `EMWindow`   | Electromagnet: Omega-reduced / A-Phi / MSC / Kelvin Benchmark | `.vol`                     |
| `radia-ih`     | `IHWindow`   | Induction Heating: 6 methods (PEEC / BEM-A / FEM A-V)         | STEP (PEEC) + `.vol` (BEM/FEM) |
| `radia-pcb`    | `PCBWindow`  | PCB / FastHenry impedance sweep                                | FastHenry `.inp`           |
| `radia-heat`   | `HeatWindow` | Thermal post-processing (q_surf -> T field on workpiece)       | `.vol` + `.sol` (q_surf)   |

All four are top-level PySide6 `QMainWindow` apps.  Launching them
spawns a single Python process (the Layer 3 panel); when the user
clicks **Run** that process spawns a Layer 4 `calc_*.py`
subprocess.  Cubit is not in either chain.

============================================================
## vol_sources — where .vol files come from
============================================================

The standalone panels accept any Netgen `.vol` (text format) that
declares the labels each method expects.  Source choices, in order
of typical lab use:

1. **Coreform Cubit** — `radia_export netgen "model.vol" order N`
   (the production source for accelerator magnets and most IH
   workpieces).  Requires Cubit licence + the radia plugin.
2. **NGSolve OCC standalone** — `OCCGeometry("model.step")` →
   `geo.GenerateMesh(maxh=...)` → `mesh.ngmesh.Save("model.vol")`.
   No Cubit licence; pure Python; meshes a STEP file from build123d
   / Onshape / FreeCAD / etc.  See `examples/` for templates.
3. **build123d → STEP → NGSolve** — the BEM-A test fixtures
   (`tests/bem/fixtures/`) use build123d to generate STEP +
   `OCCGeometry(...).GenerateMesh(...)` to produce `.vol`.
4. **Existing `.vol` from a colleague** — Kubota receives a `.vol`
   over the W: NAS share or by email; just pass it to the launcher.

The `.vol` is text and can be inspected with any editor; the
required label sets are documented per-method in the panel's
`check_method_requirements()` (e.g. for FEM A-V the `.vol` must
have `coil` / `source` / `sink` / `sibc` / `kelvin` labels).

============================================================
## vs_cubit — standalone launcher vs Cubit Solve menu
============================================================

The Cubit launcher (`Solve -> Radia-NGSolve` inside Coreform Cubit)
and the standalone .exe entry-points run the **same Layer 3 panel
window** and the **same Layer 4 calc_*.py subprocess**.  The only
real difference is what is exposed BEFORE the panel opens:

|                                            | Cubit Solve menu | Standalone `.exe` |
|--------------------------------------------|------------------|-------------------|
| Run a .jou file → mesh → `.vol` in one go | Yes              | No                |
| Run an analysis on an existing `.vol`     | Yes              | Yes               |
| Cubit licence required                     | Yes              | No                |
| Mesh inspection in 3D                      | Yes (Cubit GUI)  | No (use GMSH)     |
| `.vol` label edit / fix                    | Yes (re-export)  | No (text edit)    |

For pure analysis ("I have the `.vol`, run the panel"), the two
launchers are equivalent.  The Cubit launcher's added value is the
mesh-generation pipeline (`.jou` → Cubit batch → `radia_export
netgen` → `.vol`).

============================================================
## ih_methods — 6 IH methods on radia-ih (with/without workpiece)
============================================================

Since 4.25.0 the IH panel offers 6 methods.  Two pairs differ only
in the coil solver (PEEC perimeter filaments vs BEM-A surface RWG):

| Method                                | Coil          | Workpiece          | Inputs                  | Cost     |
|---------------------------------------|---------------|--------------------|-------------------------|----------|
| PEEC inductance (vacuum)              | PEEC filament | none               | STEP only               | ~30 s    |
| BEM-A inductance (vacuum)             | BEM-A surface | none               | STEP only               | ~25 s    |
| PEEC + BEM weak coupling              | PEEC filament | scalar BEM-SIBC    | STEP + `.vol`           | ~3 min   |
| BEM-A + BEM weak coupling             | BEM-A surface | scalar BEM-SIBC    | STEP + `.vol`           | ~3-5 min |
| PEEC coil + FEM wp (SIBC) + Kelvin    | PEEC filament | volumetric FEM SIBC | STEP + `.vol` (Kelvin) | ~4-8 min |
| Full simulation (FEM A-V + wp SIBC)   | volumetric    | volumetric FEM SIBC | `.vol` (full labels)    | ~1-7 min |

Vacuum modes (top 2) need only a coil STEP — **no `.vol` required**.
This is the cheapest way to get a coil L estimate: open the panel,
Browse the STEP, click Run, get L_coil + R_coil in <30 s.

For weak-coupled and FEM modes, the `.vol` is required and must
have the right labels (`sibc` boundary for BEM-SIBC; `coil` /
`source` / `sink` / `sibc` / `kelvin` for FEM A-V full).

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
