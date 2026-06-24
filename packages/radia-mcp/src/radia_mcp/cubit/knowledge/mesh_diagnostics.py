"""
Mesh diagnostics / cleanup / quality knowledge base.

Covers the foundational Cubit workflow for diagnosing and fixing geometry
and mesh problems. This is the "why" behind the mechanics that
cubit_scripting_knowledge shows as recipes.

Source: Coreform webinar "Your first 15 minutes with Coreform Cubit"
        (2025.8 release). Augmented with LAB field experience from the
        Radia + NGSolve workflow.
"""

DIAGNOSTICS_OVERVIEW = """
# Mesh Diagnostics & Cleanup Overview

Cubit is "semi-automated": it can recognize how to mesh a given geometry
*once the geometry is clean*. The bulk of practical meshing time is spent
in the prepare-geometry loop, not in scheme/size tuning.

The canonical beginner flow (Coreform "first 15 minutes" webinar):

```
Import CAD (STEP / native)
        |
        v
[Power Tools -> Item Wizard]
        |
        v
Check diagnostics --+-- Remove small features (blend chain, small surface)
        |           +-- Find overlaps & gaps
        |           +-- Force ACIS tolerance
        v
Preview mesh (hex hint or tet hint)
        |
        v
Mesh + validate quality (shape >= 0.2)
        |
        v
Blocks / sidesets / nodesets
        |
        v
Export
```

## When to use Power Tools vs. Command Panel

- **Power Tools + Item Wizard**: great for beginners or new geometries.
  It walks you through import -> diagnose -> clean -> mesh -> BC -> export
  and highlights the next button to click at each step.
- **Command Panel**: the fully-featured equivalent once you know what
  you want. Every Power Tool action has a Command Panel equivalent.
- **Journal + Python scripting**: reproducibility. Use the History tab
  to capture commands, toggle the Python icon to convert to
  `cubit.cmd(...)` calls, then edit and replay.

## Bottomless well: "just fix it at the CAD level"

Coreform's official recommendation is that non-watertight models or
severe overlaps should be fixed with the design engineer upstream in
the native CAD tool, not with Cubit's bottom-up repair. Cubit has
stitch / join-vertices-to-curve / force-watertight tools, but they are
harder to get right than re-exporting a clean STEP.
"""

IMPRINT_MERGE_THEORY = """
# Imprint + Merge: Why They Are Necessary

This is the single most important concept for building contiguous
(conformal) multi-body meshes in Cubit. Without it, two bodies that
share a face will mesh with *duplicate* nodes along that face, and the
resulting mesh is not continuous.

## The webcut / multi-body problem

Suppose you have two volumes sharing a face. Cubit treats each volume's
face as a *separate* surface - they are two independent topological
entities that happen to occupy the same 3D region. When you mesh:

- Volume A meshes its face with its own set of nodes.
- Volume B meshes its face with a *different* set of nodes.
- The two node sets are not linked; the solver sees a crack.

## Imprint: copy topology across volumes

`imprint volume all` projects every curve on every volume onto every
other volume. New curves and surfaces appear wherever the projection
intersects. After imprint, the *topology* matches on both sides of the
shared face - each has the same curves and sub-surfaces.

But the two matching surfaces are still separate topological entities.
Meshing them produces geometrically coincident nodes that are logically
distinct.

## Merge: identify coincident topology

`merge all` walks the topology and says "these two surfaces occupy the
same 3D region - collapse them to a single surface." After merge, there
is exactly one surface, and the mesh nodes on that surface are shared
by both volumes.

**After merge, the curve that was the shared edge disappears from the
graphics (draw). That is correct - there is only one curve, not two.**

## The full conformal pattern

```python
cubit.cmd("webcut volume 1 sweep ...")          # splits into two
cubit.cmd("imprint volume all")                 # match topology
cubit.cmd("merge all")                          # stitch
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
```

## Merge tolerance - the hidden knob

If the two surfaces are geometrically coincident but with a small gap
(typical for imports from different CAD tools), plain `merge all` will
miss them. Set the tolerance explicitly:

```python
cubit.cmd("merge tolerance 1e-4")     # absolute distance
cubit.cmd("imprint all")
cubit.cmd("merge all")
```

The right tolerance is roughly half the shortest edge you care about.
Start from the tire-toolbar heuristic:

```python
ids = cubit.parse_cubit_list("curve", "all")
shortest = min(cubit.curve(cid).length() for cid in ids)
cubit.cmd(f"merge tolerance {shortest / 2}")
```

## Verifying merge: two visualization techniques

### Curve veance tool

Toggle the "curve veance" icon in the toolbar:

- **White curves** = merged (shared between two entities).
- **Blue curves** = unmerged (boundary or genuinely unshared).

External boundaries of the model should be blue; interior edges between
merged volumes should be white.

### `draw ... with is_merged`

Run a Cubit command to visualize only merged entities:

```
draw surface with is_merged=true    # show only merged surfaces
draw curve with is_merged=false     # only unmerged curves
```

This is invaluable on complex assemblies where curve veance is visually
cluttered.

### Python API query

```python
# Find every unmerged surface (likely boundary)
sids = cubit.parse_cubit_list("surface", "all with is_merged=false")

# Find curves that SHOULD be merged but are not
unmerged_curves = cubit.parse_cubit_list(
    "curve", "all with is_merged=false"
)
```

## Pitfalls

1. **Moving a merged volume** breaks the merge. Use `move volume N z 1
   include_merged` to move both sides together.
2. **`imprint` after `mesh`** forces a full remesh - cheap commands
   should run before expensive ones.
3. **Over-merging** can happen if tolerance is too loose (merges
   unrelated surfaces). Check with `draw surface with is_merged` after
   each merge.
4. **Virtual geometry** (composite surfaces, collapsed edges) is not
   mergeable - regenerate geometry without composites first if possible.

## Related tools in other MCP modules

- `cubit_docs topic=scripting_boolean_policy` - LAB guidance on
  imprint+merge vs. unite (avoid unite for loft chains).
- `cubit_scripting_knowledge` topic `acis_tolerance` -
  ACIS-level gotchas that look like merge failures but aren't.
"""

EXODUS_DATA_MODEL = """
# Exodus Data Model: Block / Sideset / Nodeset

Cubit is built on the Exodus II mesh format (developed by Sandia, open
source, on GitHub as part of the CCAs project). Understanding how
blocks, sidesets, and nodesets differ is essential when exporting to
Abaqus, Ansys, Nastran, OpenFOAM, or the LAB Radia + NGSolve stack.

## Quick reference

| Entity   | Contains               | Preserves order?| Typical use               |
|----------|------------------------|-----------------|---------------------------|
| Block    | Elements (via volume / surface / curve or explicit) | No | Material / element formulation assignment |
| Sideset  | Element faces or edges with orientation | **Yes**  | Pressure BC, flux, traction, integrals with normal direction |
| Nodeset  | Bare nodes             | No              | Displacement / temperature / Dirichlet BC |

## Block: materials and element formulation

A block is the coarsest grouping - a set of elements that share a
material or element formulation. In Abaqus, blocks map to `*SOLID
SECTION`. In Ansys, they map to MAT / TYPE. In Radia's mesh export,
blocks are the ONLY entity translated to the `.vol` / `.msh` file.

```python
cubit.cmd("block 1 add volume 1 2 3")
cubit.cmd('block 1 name "steel_core"')
cubit.cmd("block 2 add volume 4")
cubit.cmd('block 2 name "copper_coil"')
```

Blocks can contain geometry (volumes/surfaces/curves) or raw elements
(`block 1 add tet all`). Both forms work; geometry-based blocks are
more robust under re-meshing.

## Sideset: order-preserving sub-entities for integrable BCs

A sideset is an ordered list of element faces (or edges). Order matters
because integrable boundary conditions - pressure, heat flux,
traction - need a consistent outward normal direction per face.

Sidesets are what you want for:

- Pressure loads (`p` integrated over surface)
- Flux boundary conditions (`q . n` requires normal direction)
- Periodic boundary conditions matching in pairs
- Surface-to-surface contact definitions

```python
cubit.cmd("sideset 1 add surface 5")
cubit.cmd('sideset 1 name "pressure_top"')
cubit.cmd("sideset 2 surface with x_coord < 0")    # programmatic
cubit.cmd('sideset 2 name "inlet"')
```

### Common trap: mixing surfaces and volumes in sidesets

```python
# WRONG - sidesets are surface/curve entities only
cubit.cmd("sideset 1 add volume 1")                # invalid

# CORRECT
cubit.cmd("sideset 1 add surface with volume 1")   # programmatic selector
```

## Nodeset: bare nodes for Dirichlet BCs

A nodeset is a set of nodes without connectivity information. Use it for:

- Displacement / velocity boundary conditions (applied at nodes)
- Temperature / voltage Dirichlet conditions
- Point loads
- Locking constraints

```python
cubit.cmd("nodeset 1 add surface 1")
cubit.cmd('nodeset 1 name "clamped_base"')
cubit.cmd("nodeset 2 add vertex 1 2 3")
cubit.cmd('nodeset 2 name "corners"')
```

### Nodeset vs. sideset: quick decision

- "I have a pressure in MPa" -> sideset (need normal direction).
- "I have a displacement in mm" -> nodeset.
- "I have a heat flux in W/m^2" -> sideset.
- "I have a fixed temperature in K" -> nodeset.

## LAB policy: `export` uses blocks only

The `export netgen` / `export gmsh` / `export jmag_nastran`
commands (from the `cubit_mesh_export` LAB plugin) read from blocks
only. Sidesets and nodesets are ignored by this plugin.

If you need boundary conditions in the downstream solver:

- For **NGSolve**: name the *OCC faces* before meshing, not sidesets.
  See `cubit_docs topic=netgen` for the `set_name_for_occ_face` workflow.
- For **Abaqus / Nastran / Gmsh v4.1**: use plain `export abaqus` etc.,
  which do translate sidesets and nodesets.

See `cubit_docs topic=scripting_blocks_only_policy` for details.

## Naming convention

Always name every block, sideset, and nodeset. Downstream solvers and
post-processors key off names, and a name collision between blocks and
sidesets (both use ID 1 for example) can silently overwrite data.

```python
# SAFE: distinct namespace prefixes
cubit.cmd('block 1 name "vol_steel"')
cubit.cmd('sideset 1 name "srf_pressure"')
cubit.cmd('nodeset 1 name "pts_clamped"')
```
"""

SMALL_FEATURE_CLEANUP = """
# Small Feature Cleanup: detect, remove, re-check

Small features (tiny fillets, thin flanges, sub-millimeter holes) are
the #1 reason Cubit cannot generate a meshable topology on imported
STEP. The cleanup loop is:

```
detect small features
   |
   v
remove (blend chain for fillets; remove surface for small patches)
   |
   v
re-detect
   |
   v
loop until no small features remain
```

## Using Power Tools (beginner-friendly)

In the Power Tools panel (`View -> Power Tools`), click the Item Wizard
hat icon. Step through to "Prepare Geometry". Click "Remove small
features". A dialog lets you set a threshold (rule of thumb: 5% of the
typical element size - if your mesh size is 1.0, threshold is 0.05).

After detection, a list appears. Click a curve / surface to highlight
it in the graphics window. Click **Remove blend chain** for fillets or
**Remove surface** for small patches. The Power Tool issues the
corresponding Command Panel command, so you can copy it to a journal
file.

## Command-panel equivalents

```python
# Detect small curves / surfaces
cubit.cmd("find small features curves size 0.05")
cubit.cmd("find small features surfaces size 0.01")

# Remove blend chain (fillet removal - propagates along the fillet chain)
cubit.cmd("remove surface 42 extend")

# Remove a small flat surface (absorbs into neighbor)
cubit.cmd("remove surface 77")

# More aggressive: remove small features globally
cubit.cmd("autoclean volume all small_details yes small_curves yes")
```

## When blend chain removal is right

- Rolling fillets along an edge that wrap a whole corner - `remove
  surface <id> extend` finds and removes the whole chain in one call.
- Do NOT blanket-remove fillets if they are load-bearing (stress
  concentrations matter); check with the design engineer first.

## When to refuse to clean

- Thin flanges critical to structural behavior - mesh them with a thin
  shell element type (`quad`/`tri` with thickness attribute) instead of
  merging away the feature.
- Small holes used for fluid inlets - use a finer local size rather
  than removing the hole.

## After cleanup: re-run check diagnostics

```python
cubit.cmd("check diagnostics")
```

The Power Tools panel surfaces this as a button. Output includes:

- Watertightness (should be OK after cleanup).
- Meshability (blocks / volumes categorized as hex-able / sweep-able /
  tet-only).
- Overlaps / gaps (see next topic).
"""

GAPS_AND_OVERLAPS = """
# Gaps & Overlaps: detection and repair

When multiple CAD bodies are assembled (e.g., bolt through hole, gasket
between plates, or nested coil geometry), numerical roundoff often
leaves gaps or overlaps at their interfaces. Cubit provides native
detection and repair.

## Detection

```python
cubit.cmd("find overlap volume all")
cubit.cmd("find gap volume all")
```

The output lists volume pairs and the magnitude of overlap/gap. In the
Command Panel, the same workflow is under Geometry -> Cleanup -> Find
Gaps and Overlaps.

## Repair: remove overlap from one volume

```python
cubit.cmd("remove overlap volume 2 from volume 1")
# Reads: "remove the part of volume 2 that overlaps volume 1" from
# volume 2's representation. Volume 2 shrinks.
```

The direction matters - always remove the overlap from the *smaller*
or less-critical volume, so you preserve the geometry of the important
one.

## Repair: close a gap

For small gaps, re-export the STEP from the native CAD tool with a
tighter tolerance if possible. Otherwise:

```python
# Option A: loose merge (works if gap <= merge tolerance)
cubit.cmd("merge tolerance 1e-3")
cubit.cmd("imprint all")
cubit.cmd("merge all")

# Option B: create a curve bridging the gap, then a bounding surface
cubit.cmd("create curve vertex 5 6")
cubit.cmd("create surface curve 101 102 103")
cubit.cmd("stitch volume 1 2")
```

## Non-watertight STL / faceted input

STLs have their own failure modes: non-watertight, flipped normals,
tiny duplicate triangles. Coreform's recommendation:

- Heal STL upstream in MeshLab / 3D printing slicer / EOS software.
- Then import into Cubit with the merge-stitch option:

```python
cubit.cmd('import stl "model.stl" feature_angle 135.0 merge stitch')
```

- For tet-only meshes on facet geometry, Cubit can meshes STL directly
  but cannot perform robust repairs like brep cleanup.

## Workflow guidance

1. Always fix watertightness upstream (CAD) if possible.
2. If stuck with a problematic STEP, cleanup in this order:
   detect -> remove small features -> re-detect -> imprint+merge ->
   find gaps / overlaps -> repair interface pairs -> final check
   diagnostics.
3. Snapshot (`cubit_checkpoint`) before each destructive repair so you
   can roll back.

## LAB-specific guidance

For the Radia + NGSolve stack, the recommended flow when assemblies
have alignment tolerance issues is:

- Export the assembly from build123d with `export_step_with_labels` so
  the faces are pre-named (see memory
  `project_joachim_glue_step_report.md`).
- Use `cubit.cmd("import step ... no_heal")` to prevent Cubit's healing
  engine from throwing away faces.
- Apply imprint+merge with a tolerance matching the CAD export
  tolerance (typically 1e-6 m for build123d).
"""

MESH_QUALITY_DEEP_DIVE = """
# Mesh Quality: Metrics, Thresholds, Convergence

Cubit's default quality metric is `shape`, a meta-metric that rolls
together size, distortion, scaled Jacobian, and a few others into a
single 0-1 score (higher is better).

## Threshold rules of thumb

| Element type | shape (pass) | scaled_jacobian (pass) | aspect_ratio (pass) |
|--------------|--------------|------------------------|---------------------|
| Tet          | >= 0.2       | >= 0.2                 | <= 5                |
| Hex          | >= 0.3       | >= 0.5                 | <= 3                |
| Wedge        | >= 0.2       | >= 0.2                 | <= 5                |
| Pyramid      | >= 0.2       | >= 0.2                 | <= 5                |

- `>= 0.2` for tet `shape` / `scaled_jacobian` is Cubit's canonical
  "pretty good" threshold.
- Anything `< 0.05` is a degenerate element and typically breaks
  solvers (NaN Jacobian, singular mass matrix).

## Query quality

```python
# Global quality report (histogram + min/max)
cubit.cmd("quality volume all shape global")

# Per-element quality via Python API
tets = cubit.parse_cubit_list("tet", "all")
for t in tets:
    q = cubit.get_quality_value("tet", t, "shape")
    if q < 0.2:
        print(f"tet {t}: shape = {q:.3f} (bad)")
```

## Visualize quality

```python
cubit.cmd("draw tet all shape")
# renders tets colored by shape quality - dark red = bad, green = good

cubit.cmd("draw tet with shape < 0.2 color red")
# isolate only the bad elements for inspection
```

## Convergence workflow (auto-refine)

See `cubit_docs topic=scripting_quality_convergence_loop` or
`cubit_toolbar_guide topic=quality_convergence` for a full pattern.
The core idea is:

1. Mesh with initial size.
2. Query min quality.
3. If >= threshold, done.
4. Else shrink size (or refine locally near bad elements) and remesh.
5. Budget the loop (MAX_ITERATIONS <= 5) to prevent runaway.

## Local refinement vs. global size shrink

- **Global size shrink**: simple, works, but element count explodes.
  Use only for small models.
- **Local refinement**: `refine tet <id> numsplit 1` is surgical.
  Identify bad elements first, then refine only them.
- **Size function**: `volume 1 size_function <sf>` adapts size by
  geometry (curvature, small feature distance). Most efficient but
  steepest learning curve.

## Integrating with the solver

Check that the solver's Jacobian warnings align with Cubit's quality
warnings. For NGSolve (via export netgen), the high-order mesh
curving can rescue borderline straight-edge tets - so a `shape` score
just below 0.2 is often survivable if the mesh is curved at order 2+.

```python
cubit.cmd(f'export netgen "{out_path}" order 2 overwrite')
```

See `netgen_workflow_guide workflow=accuracy` for how to verify order
selection against analytical reference problems.

## Coreform Cubit 2026.6 quality note

The 2026.6 release highlights more robust triangle/tet meshing, higher
precision normalized mesh-quality metrics, and Jacobian/scaled-Jacobian
metrics for Tetra10/Tri6. Keep Cubit as the mesher's own quality oracle,
but do not make it the only gate: after `export netgen`, run Radia's
solver-neutral `.vol` intake checks as an independent confirmation.

Recommended public checks:

```powershell
python examples/cubit_mesh_export/validation_vol_tet_quality.py --vol C:\\temp\\model.vol
python examples/cubit_mesh_export/validation_vol_surface_triangle_quality.py --vol C:\\temp\\model.vol
```

`NetgenTriTetVolMesh.surface_triangle_quality_summary()` reports boundary
triangle area, edge-ratio, angle range, and `2 * inradius / circumradius`
quality. Use it as the FEM/BEM trace gate: bad surface triangles can break
scalar-BEM/RWG behavior even when the volume tet inventory looks acceptable.
"""

JOURNAL_TO_PYTHON_WORKFLOW = """
# Journal -> Python Conversion Workflow

Cubit has two scripting layers:

- **Journal files** (`.jou`): Cubit Command Language. Line-by-line
  commands. Played back with `cubit.cmd("play 'file.jou'")` or from
  the GUI File -> Play Journal.
- **Python scripts**: `cubit.cmd("...")` strings, but wrapped in Python
  so you can use f-strings, numpy, loops, try/except, and any library.

The natural workflow is:

1. Explore interactively in the GUI (click through Power Tools /
   Command Panel).
2. Capture the commands from the History tab.
3. Convert to a journal file (saves as .jou).
4. Toggle the Python icon in the journal editor - Cubit converts each
   `cmd` line to `cubit.cmd("cmd")` in Python.
5. Edit the Python script - add loops, parameters, error handling.
6. Replay as a repeatable pipeline.

## Capturing commands from History

In the Cubit GUI, the bottom-right pane has tabs:

- **Command Panel**: GUI-driven command construction.
- **Power Tools**: beginner wizard.
- **History**: log of every command issued this session.
- **Errors**: warnings / errors only.
- **Console**: interactive Python prompt inside Cubit.

Right-click the History tab -> **Save to journal file** to capture
exactly what you did.

## Toggling Python mode

With a journal editor open, there is a small Python icon in the
toolbar. Click it to switch the editor between:

- Cubit command language mode (plain commands, one per line).
- Python mode (`cubit.cmd("command")` wrapped).

The conversion is bidirectional and lossless *for simple commands*. If
you edit the Python version and add `import numpy`, you cannot toggle
back to command-language mode without losing the import.

## Canonical Python structure

```python
#!python
\"\"\"Replay of interactive session, converted from journal.\"\"\"

import os
import numpy as np

# When run outside Cubit (CI, subprocess), initialize:
if __name__ == '__main__':
    import cubit
    cubit.init(['cubit', '-nojournal', '-batch'])

cubit.cmd("reset")

# Parameterize what was hardcoded in the journal
SIZE = float(os.environ.get("MESH_SIZE", "0.1"))
ORDER = int(os.environ.get("ORDER", "2"))

# Captured commands - now parameterized
cubit.cmd("create brick x 10")
cubit.cmd(f"volume all size {SIZE}")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Block registration always required for export
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "domain"')

# Export
out = os.environ.get("OUT_PATH", "out.vol")
cubit.cmd(f'export netgen "{out}" order {ORDER} overwrite')
```

## Running scripts

| Context              | Command                                        |
|----------------------|------------------------------------------------|
| From Cubit GUI       | `cubit.cmd('play filename.py')` in console     |
| Cubit batch          | `coreform_cubit -batch -nographics -nojournal file.py` |
| Python subprocess    | `python file.py` with `cubit.init(...)`        |

## Why journal + Python instead of just Python

The journal file is a stable artifact that can be replayed inside a
Cubit session to recreate exact state. Python scripts that call `cubit.
cmd` are more flexible but embed Cubit initialization and session
state. Keep both in a project:

- `journals/build_geometry.jou` - pure Cubit commands, replayable.
- `scripts/mesh_with_refinement.py` - Python orchestration that loads
  the journal + adds logic.

See `cubit_docs topic=scripting_aprepro_journal` for APREPRO
expressions (Cubit's built-in templating) as a middle ground between
plain journal and full Python.
"""

# ---- Topic registry ----------------------------------------------------

def get_diagnostics_documentation(topic: str = "overview") -> str:
    """Return mesh diagnostics / cleanup documentation by topic."""
    topics = {
        "overview": DIAGNOSTICS_OVERVIEW,
        "workflow": DIAGNOSTICS_OVERVIEW,                   # alias
        "power_tools": DIAGNOSTICS_OVERVIEW,                # alias
        "item_wizard": DIAGNOSTICS_OVERVIEW,                # alias

        "imprint_merge": IMPRINT_MERGE_THEORY,
        "contiguous_mesh": IMPRINT_MERGE_THEORY,            # alias
        "merge_theory": IMPRINT_MERGE_THEORY,               # alias
        "curve_veance": IMPRINT_MERGE_THEORY,               # alias
        "is_merged": IMPRINT_MERGE_THEORY,                  # alias

        "exodus": EXODUS_DATA_MODEL,
        "block_sideset_nodeset": EXODUS_DATA_MODEL,         # alias
        "data_model": EXODUS_DATA_MODEL,                    # alias

        "small_features": SMALL_FEATURE_CLEANUP,
        "blend_chain": SMALL_FEATURE_CLEANUP,               # alias
        "cleanup": SMALL_FEATURE_CLEANUP,                   # alias
        "fillet_removal": SMALL_FEATURE_CLEANUP,            # alias

        "gaps_overlaps": GAPS_AND_OVERLAPS,
        "overlaps": GAPS_AND_OVERLAPS,                      # alias
        "find_overlap": GAPS_AND_OVERLAPS,                  # alias
        "stitch": GAPS_AND_OVERLAPS,                        # alias

        "quality": MESH_QUALITY_DEEP_DIVE,
        "quality_deep": MESH_QUALITY_DEEP_DIVE,             # alias
        "shape_metric": MESH_QUALITY_DEEP_DIVE,             # alias
        "scaled_jacobian": MESH_QUALITY_DEEP_DIVE,          # alias

        "journal_to_python": JOURNAL_TO_PYTHON_WORKFLOW,
        "history_tab": JOURNAL_TO_PYTHON_WORKFLOW,          # alias
        "journal": JOURNAL_TO_PYTHON_WORKFLOW,              # alias
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join([
            DIAGNOSTICS_OVERVIEW,
            IMPRINT_MERGE_THEORY,
            EXODUS_DATA_MODEL,
            SMALL_FEATURE_CLEANUP,
            GAPS_AND_OVERLAPS,
            MESH_QUALITY_DEEP_DIVE,
            JOURNAL_TO_PYTHON_WORKFLOW,
        ])
    if topic in topics:
        return topics[topic]
    return (
        f"Unknown topic: '{topic}'. "
        f"Available: all, overview, imprint_merge, exodus, "
        f"small_features, gaps_overlaps, quality, journal_to_python. "
        f"(Aliases accepted: curve_veance, is_merged, "
        f"block_sideset_nodeset, blend_chain, find_overlap, "
        f"shape_metric, history_tab.)"
    )
