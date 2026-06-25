# build123d → Netgen → Gmsh pipeline scaffold

Minimal, loopable scaffold exercising the lab's main CAD/mesh/post stack:

```
build123d (OCCT)  →  Netgen (tet via netgen.occ)  →  Gmsh .msh/.geo files
    CAD                    mesh                              visualization
```

**Purpose**:
1. Keep a known-good end-to-end integration under version control.
2. Serve as **training input for the lab mcp-servers**
   (`mcp-server-build123d`, `mcp-server-gmsh`) — every run produces a
   `.brep / .msh / _post.msh / .json` quartet that Claude Code can inspect
   and reason about.
3. Be a drop-in skeleton that research projects can fork.

## Lab policy reminder (2026-04-19)

- **CAD**: build123d (Python, OCCT). FreeCAD is not used on dev machines.
- **Mesh (tet, main)**: Netgen via `netgen.occ.OCCGeometry`.
- **Mesh (hex, sub)**: Cubit (required for Radia/ELF hex). Not exercised here.
- **Post**: write GMSH `.msh` data blocks directly (`$NodeData`,
  `$ElementData`, `$PhysicalNames`) and use standalone GMSH only for display.
  The public examples do not import the pip `gmsh` runtime.
- **Mesh generation through Gmsh is forbidden** (covered by the
  `mcp-server-gmsh` lint rules).

See `docs/research/policy/strategy.md` and `toolchain.md` for the full rationale.

## Files

| File | Purpose |
|---|---|
| `_pipeline.py` | Core helpers: `run_pipeline(part, …)` single-region, `run_pipeline_multi(regions, …)` multi-region |
| `demo_box.py` | Minimal one-case demo (plate with hole, single region) |
| `demo_ih_multi.py` | Multi-region demo (workpiece + coil + air) |
| `sweep.py` | 19-case parametric sweep across 4 geometry families (3 single + 1 multi) |
| `validation_helix_mesh_sweep.py` | Heavier validation-class helix conductor sweep; records analytic/CAD volume error plus mesh/post stats |
| `validation_halbach_region_sweep.py` | Heavier validation-class segmented Halbach sweep; checks CAD volume, magnetization labels, and one mesh region per segment |
| `validation_coaxial_region_stack.py` | Validation-class touching coaxial regions; checks analytic shell volumes and region-name preservation through STEP/Netgen/Gmsh |
| `validation_laminated_stack_region_sweep.py` | Validation-class touching laminated box stack; checks per-layer volumes, fill factor, and named region preservation |
| `validation_racetrack_plate_air_region.py` | Validation-class racetrack coil + conductive plate + air box; checks analytic region volumes and named region preservation |
| `validation_build123d_cubit_measurement.py` | Validation-class build123d STEP round-trip measured by headless Cubit API; checks volume, surface-area, and bbox parity |
| `validation_build123d_measurement_health.py` | Validation-class build123d assembly measurement health; reports volume fractions, bbox fill fraction, and worst Cubit volume/area/bbox mismatches |
| `validation_build123d_parameter_sweep_summary.py` | Validation-class build123d parameter sweep design table; checks monotonic volume/area trends and constraint violations before meshing |
| `validation_build123d_bbox_clearance_audit.py` | Validation-class build123d bbox clearance audit; separates provably disjoint pairs from all-axis bbox overlap pairs that need precise geometry checks |
| `validation_enclosure_cubit_measurement.py` | Validation-class enclosing-box/void-region STEP round-trip measured by headless Cubit API; checks bbox margin, analytic volume/area, and Cubit volume/area/bbox parity |
| `validation_build123d_cubit_boundary_normals.py` | Validation-class box boundary normals; checks analytic build123d face vector areas, optionally against a named Netgen `.vol` boundary mesh |
| `validation_build123d_cubit_pressure_force.py` | Validation-class box pressure force; checks analytic build123d face forces, optionally against named Netgen `.vol` pressure-force rows |
| `validation_build123d_cubit_pressure_moment.py` | Validation-class box pressure force/moment; checks analytic build123d face moments, optionally against named Netgen `.vol` pressure-moment rows |
| `validation_build123d_cubit_pressure_resultant.py` | Validation-class box pressure resultant; checks analytic build123d face force/moment summary against named Netgen `.vol` boundary summary |
| `validation_build123d_cubit_traction_moment.py` | Validation-class box vector-traction force/moment; checks analytic build123d face moments, optionally against named Netgen `.vol` vector-traction rows |
| `runs/` | Output directory (`*.brep` / `*.step`, `*.msh`, `*_post.msh`, `*.json`, `sweep_summary.json`) |

## Run

```
cd examples/build123d_netgen_gmsh_flow
python demo_box.py              # single-region demo
python demo_ih_multi.py         # multi-region demo (workpiece/coil/air)
python sweep.py --quick         # 4 cases (one per family, smoke)
python sweep.py                 # full sweep (19 cases)
python validation_helix_mesh_sweep.py --quick  # validation-class baseline
python validation_helix_mesh_sweep.py          # heavier 4-case helix sweep
python validation_halbach_region_sweep.py --quick
python validation_halbach_region_sweep.py
python validation_coaxial_region_stack.py --quick
python validation_coaxial_region_stack.py
python validation_laminated_stack_region_sweep.py --quick
python validation_laminated_stack_region_sweep.py
python validation_racetrack_plate_air_region.py --quick
python validation_racetrack_plate_air_region.py
python validation_build123d_cubit_measurement.py --require-cubit
python validation_build123d_measurement_health.py --require-cubit
python validation_build123d_parameter_sweep_summary.py
python validation_build123d_bbox_clearance_audit.py
python validation_enclosure_cubit_measurement.py --require-cubit
python validation_build123d_cubit_boundary_normals.py
python validation_build123d_cubit_boundary_normals.py --vol C:\temp\box.vol --out C:\temp\box_boundary_normals_summary.json
python validation_build123d_cubit_pressure_force.py
python validation_build123d_cubit_pressure_force.py --vol C:\temp\box.vol --out C:\temp\box_pressure_force_summary.json
python validation_build123d_cubit_pressure_moment.py
python validation_build123d_cubit_pressure_moment.py --vol C:\temp\box.vol --out C:\temp\box_pressure_moment_summary.json
python validation_build123d_cubit_pressure_resultant.py
python validation_build123d_cubit_pressure_resultant.py --vol C:\temp\box.vol --out C:\temp\box_pressure_resultant_summary.json
python validation_build123d_cubit_traction_moment.py
python validation_build123d_cubit_traction_moment.py --vol C:\temp\box.vol --out C:\temp\box_traction_moment_summary.json
```

On a warm Python (all imports cached) the full sweep takes tens of seconds.
The helix validation is deliberately not a pytest test; it is a heavier
example that writes reusable records under `runs/validation_helix_mesh_sweep/`.

## Pipeline stages

### 1. CAD — build123d → BREP
- `export_brep(part, "<label>.brep")`
- record: `volume`, `area`, `faces`, `edges`, `min_edge`, `bbox_size`,
  `is_valid`
- units: build123d default (mm)

### 2. Mesh — Netgen tet
- `OCCGeometry("<label>.brep")` → `.GenerateMesh(maxh=...)`
- export via `NgMesh.Export("<label>.msh", "Gmsh2 Format")`
  (Gmsh v2.2 ASCII; v4 is also available as `"Gmsh Format"` if needed)
- record: `nv / ne / nface / nedge`, `gen_seconds`, `maxh`

### 3. Post — GMSH file writer
- reads Netgen's exported node and element tags from `<label>.msh`
- builds a **dummy scalar field** `f(x,y,z) = x + 2y + 3z`
  (placeholder for real solver output)
- writes `<label>_post.msh` as GMSH v2.2 ASCII with a `$NodeData` block
- does not import the pip `gmsh` runtime; standalone GMSH can open the file
  later for display

## Multi-region flow (`run_pipeline_multi`)

**Why it's tricky**: build123d's STEP export does not emit region names in
a form Netgen picks up, and Netgen's Gmsh-v4 exporter preserves per-solid
physical *tags* but not *names*. We bridge the three-step gap explicitly:

1. **build123d side**: pack regions into a `Compound(children=[...])` and
   export as STEP. Children order is the contract — it becomes the region
   index.
2. **Netgen side**: load the STEP, iterate `geo.shape.solids` in order,
   call `s.mat(name)`, then `Glue(sols)` and rebuild the `OCCGeometry`.
   `Glue` is what makes Netgen treat the named solids as separate domains
   that survive to meshing. `ngsolve.Mesh(ng_mesh).GetMaterials()` now
   returns the right names.
3. **GMSH file side**: read Netgen's exported `.msh`, sort volume physical
   tags, write `$PhysicalNames`, and emit an `$ElementData` view called
   `<label>_region_id` where each tet carries its 1-based region index.
   Standalone GMSH renders each region a different color by default.

Usage:

```python
from _pipeline import run_pipeline_multi

regions = [
    (workpiece_part, "workpiece"),
    (coil_part,      "coil"),
    (air_part,       "air"),
]
rec = run_pipeline_multi(regions, out_dir="./runs", label="ih_multi",
                         maxh=6.0)
```

Output extras (multi mode):
- `<label>.step` — Compound with all regions
- `<label>.msh` — Netgen Gmsh-v4 output (physical tags only)
- `<label>_post.msh` — **named** physical groups + region_id view
- `<label>.json` — per-region element counts, physical tags, names

Region-order is a load-bearing contract — the list you pass to
`run_pipeline_multi` must match the spatial containment you want in the
mesh (usually: small inner regions first, large surrounding region last).

## Extending

### Add a geometry family to the sweep

1. Add a builder function in `sweep.py` returning a build123d Shape
   (set `.label`).
2. Add a generator `XXX_grid()` yielding dicts with `builder`, `kwargs`,
   `label`, `maxh`.
3. Append the generator to `all_cases()`.

### Replace the dummy field with real solver data

In `_pipeline.py::_stage_post`, swap the `values = [...]` expression for
data from your solver (Radia B-field, NGSolve potential, etc.). Keep the
tag / data order consistent with `node_tags` from `getNodes()`.

For model-based data referencing mesh entity tags directly (rather than
node coords), write an `$ElementData` block keyed by the element tags read
from the exported `.msh`.

### Replay / inspect via mcp-server

After a run, the `.json` record is a compact summary; the `.brep` +
`_post.msh` are browsable:

- `mcp-server-build123d`: `inspect_geometry("runs/<label>.brep")` for CAD
  quality warnings (micro-edges, non-valid shapes, face/edge histograms).
- `mcp-server-gmsh`: `gmsh_usage("workflow")` and `gmsh_reference("msh_format")`
  for display-file conventions and the GMSH data-block format.

## Known limitations / TODO

- Single-region pipeline writes `.msh` in Gmsh v2.2 (2nd order cap).
  Multi-region pipeline writes v4 ("Gmsh Format") because v2.2 flattens
  physical groups.  Bumping the single-region path to v4 is cheap if
  higher-order elements are needed later.
- **Region contract is positional, not by label**. build123d-side labels
  are collected but don't propagate through STEP; the `(part, name)`
  tuple order is the single source of truth. If you add/reorder regions,
  re-check every consumer.
- `run_pipeline_multi` assumes **boolean-disjoint regions** (typically
  built as `outer - inner1 - inner2`). Overlapping regions confuse
  `Glue`. Always subtract inner regions from the outer air domain.
  `shape_envelope_row`, `enclosing_box`, `enclosure_clearance_row`, and
  `enclosure_difference_region` make that contract explicit before export.
- No hex path. Cubit (`.jou`) path is a separate scaffold (not yet built).
- Real solver data (Radia B-field, NGSolve potential) is not wired up —
  see the "Replace the dummy field" note above.
