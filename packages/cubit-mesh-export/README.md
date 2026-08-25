# cubit-mesh-export

Solver-neutral mesh export from [Coreform Cubit](https://coreform.com/products/coreform-cubit/) to [NGSolve](https://ngsolve.org)/[Netgen](https://github.com/NGSolve/netgen).

`cubit-mesh-export` is the **shared infrastructure layer** in the Radia
toolchain. It ships mesh export, the Kelvin open-boundary
transformation, symmetry helpers, and the Dirichlet label conventions
that every domain-specific Radia notebook or headless workflow consumes.

## Features

- **Cubit plugin** (`.ccm` + `.pyd`, Coreform Cubit 2025.12+):
  - `export {netgen|gmsh|vtk|femeem|meg|nastran_bdf}` APREPRO commands
  - **Export Mesh** GUI menu / toolbar inside Cubit's embedded Python
- **Arbitrary-order curving** (order 1-5) via ACIS geometry projection
- **Kelvin open-boundary** transformation built into `export netgen`
  (auto-add an exterior sphere with copy-mesh + periodic identification)
- **Per-axis symmetry-plane BC labels** (`bn`/`ht`) for 1/2 and 1/4 reduced
  domains
- **Dirichlet / Neumann label conventions** at three levels
  (BND / BBND / BBBND -- see table below)
- **Companion JSON** beside every `.vol` with CAD reference values for
  Volume / Area / Length consistency checking
- **Free-mesh sideset preservation**: direct tri/quad faces from Sculpt and
  imported Exodus sidesets are exported even when they have no parent CAD
  surface, with duplicate geometry-owned faces removed
- **Standalone checker** that does NOT require Cubit (`check-vol` CLI)

## Install

```bash
pip install "radia[cubit]"
cubit-plugin-install
```

The second command deploys the Cubit plugin binaries, the Netgen DLLs,
the Cubit-side Python helpers (`cubit_helpers/add_kelvin.py`,
`cubit_helpers/auto_kelvin_entry.py`), and the Radia Export Mesh toolbar
startup registration into your Coreform Cubit 2025.12 profiles.  The toolbar
runs only inside Cubit's embedded Python; normal Radia Python uses notebooks
and headless scripts and does not need PySide6.
Use `cubit-plugin-install --all-users` for a shared lab machine.

### Upgrade

```bash
pip install --upgrade cubit-mesh-export
cubit-plugin-install
```

Always re-run `cubit-plugin-install` after upgrading.
`cubit-plugin-install --verify-only` checks both the deployed binary
hashes and, when `radia` is installed, the Cubit toolbar startup
registration.

## Cubit commands

```
export netgen "model.vol" order 3 overwrite                 # NGSolve FEM (.vol)
export gmsh   "model.msh" order 2 overwrite                 # GMSH v4.1 raw data + .geo launch
export nastran_bdf "model.bdf" order 2 overwrite             # Nastran BDF
export vtk    "model.vtk" order 2 overwrite                 # VTK Legacy
```

`export jmag_nastran` remains as a deprecated compatibility alias for old
journals.  Nastran output is a mesh-interchange deck: blocks become PSOLID or
PSHELL properties, sidesets become collision-free PSHELL properties, and
nodesets become SET1 cards.  The exporter does not invent MAT cards; assign
physical material data in the receiving solver.

For Radia post-processing, `.geo` is the standard launch artifact. The Gmsh
export writes:

- `model.msh`: raw GMSH v4.1 mesh/data container
- `model.geo`: normal review entry point; it merges `model.msh`
- `model.geo.opt`: exact Gmsh sidecar auto-loaded when `model.geo` opens
- `model.msh.opt`: raw mesh/data inspection sidecar when `model.msh` opens

Associate/open `.geo` for normal review; treat `.msh` as optional raw
mesh/data inspection.

The `export netgen` command additionally accepts Kelvin / symmetry
options (see below). The other formats do not consume Kelvin.

### Free Sculpt meshes

Sculpt places its hexes and boundary faces in free mesh groups.  Use
`gen_sidesets 2` and put the hexes in a material block before export:

```text
sculpt volume all processors 1 size 0.01 gen_sidesets 2
block 1 add hex all
block 1 name "iron"
export netgen "model.vol" order 1 overwrite
```

The Netgen exporter reads both geometry-owned surfaces and direct/free
sideset tri/quads, deduplicates equal connectivity, and emits a synthetic
boundary descriptor for the latter.  It recovers `DomainIn`/`DomainOut` from
the adjacent block elements and splits one sideset into multiple descriptors
when it spans more than one material pair.  This is essential for disconnected
free bodies and mixed-material interfaces: assigning every free face to domain
1 silently removes the other materials from boundary-based operators.
`check-vol` audits that every volume domain appears in a valid boundary or
interface descriptor.  Also verify that the `.vol` surface-element count equals
an independent topological-skin count before a surface-charge formulation uses
the mesh.  Sculpt material blocks can exist without an owning CAD volume.  The
companion JSON records these under `mesh_only_materials`, rather than inventing
a zero CAD volume in `materials`; `check-vol` then treats the label as valid but
reports that no CAD-volume reference is available.

For a free-mesh material interface created after Sculpt, materialize the skin
of one block as a sideset:

```text
skin block 1 make sideset 3
sideset 3 name "left_right_interface"
```

The exporter removes connectivity already owned by the exterior sideset,
emits the remaining interface once with both nonzero domains, and computes the
sidecar area from the faces actually exported under that label.  `check-vol`
reports exterior and internal-interface counts separately and rejects duplicate
surface connectivity.

## Workflow

```
   ┌────────────┐    export netgen        ┌──────────┐
   │  Cubit     │   ─────────────────────────▶  │  .vol    │
   │  geometry  │   (+add_kelvin, +sym)         │          │
   └────────────┘                                └──────────┘
                                                      │
                                                      ▼
                       user opens a Radia Simulink application block;
                       IH may also use its temporary comparison notebook
```

`cubit-mesh-export` produces the `.vol` and the label conventions; the
domain block/headless tool reads the `.vol` and applies the physics. There is no
"pick your analysis" launcher in this plugin -- end-user tools split by
analysis target (IH designer / electromagnet designer / ...), not by
solver type.

## Kelvin open-boundary transformation

Idempotent helper: skipped if a `kelvin` block already exists; needs an
`air` block in the current Cubit model.

```
export netgen "model.vol" order 3 overwrite \
    add_kelvin                       # auto-create the exterior Kelvin sphere
    [kelvin_air "air"]               # name of the air block (default "air")
    [kelvin_block "kelvin"]          # name to give the Kelvin block (default "kelvin")
    [kelvin_mesh 0.03]               # tet size [m] on the Kelvin shell
                                     # (omit to inherit from air outer surface)
    [kelvin_sym_x {off|bn|ht}]       # per-axis symmetry-plane BC
    [kelvin_sym_y {off|bn|ht}]       # off = no reduction (default)
    [kelvin_sym_z {off|bn|ht}]       # bn  = B.n=0    (flux parallel)
                                     # ht  = HxN=0    (flux perpendicular)
```

In the Export Mesh GUI, the same options appear as widgets on the
Netgen Vol export dialog (only there -- Kelvin is `.vol`-specific).

The Kelvin step runs **before** the mesh extract / .vol write, so the
new `kelvin` block, the `kelvin_int` / `kelvin_ext` sidesets, and the
optional `sym_<bc>_<axis>` sidesets all end up in the `.vol`.

### Symmetry semantics

| `kelvin_sym_<axis>` | Sideset name produced  | B/H constraint    | Radia image | A formulation | Omega formulation |
|--------------------|------------------------|-------------------|-------------|---------------|-------------------|
| `off`              | (none)                 | (full domain)     | n/a         | n/a           | n/a               |
| `bn`               | `sym_bn=0_<axis>`      | B·n = 0           | `+`         | Dirichlet (A×n=0) | natural          |
| `ht`               | `sym_ht=0_<axis>`      | H×n = 0           | `-`         | natural          | Dirichlet (Ω=const) |

The convention is **physics-named, formulation-agnostic**: the same
`sym_bn=0_x` sideset means "B.n = 0 on x = 0 plane" regardless of
whether the domain panel solves A or Omega. Each domain tool decides
which BC type to apply per its formulation.

1/8 reduction (all three axes set to `bn` or `ht`) is supported when at
least one axis is `ht`. Three `bn` axes are physically impossible (B
parallel to three mutually perpendicular planes forces B = 0
everywhere) and rejected.

## Label conventions

`cubit-mesh-export` reserves a small set of label names and prefixes
across all three NGSolve dimension levels (BND / BBND / BBBND). The
domain tools rely on these to wire up Dirichlet / Kelvin / symmetry
without having to inspect geometry.

### BND -- surface labels (NGSolve `mesh.GetBoundaries()`)

Source: Cubit **sidesets** on surfaces.

| Cubit sideset name  | NGSolve BND name    | Meaning                                                     |
|---------------------|---------------------|-------------------------------------------------------------|
| `kelvin_int`        | `kelvin_int`        | Inner Kelvin face (auto-paired with outer via copy-mesh)    |
| `kelvin_ext`        | `kelvin_ext`        | Outer Kelvin face                                           |
| `sym_bn=0_<axis>`   | `sym_bn=0_<axis>`   | B.n = 0 (flux parallel) symmetry plane                      |
| `sym_ht=0_<axis>`   | `sym_ht=0_<axis>`   | H×n = 0 (flux perpendicular) symmetry plane                 |
| `dir_<name>`        | `dir_<name>`        | **Dirichlet** surface (variable = 0; physics is solver-side)|
| `neu_<name>`        | `neu_<name>`        | Explicit Neumann (= no-op; documentation only)              |
| anything else       | (passes through)    | Free-form name; meaning is up to the domain tool            |

`kelvin_int` / `kelvin_ext` are auto-detected from the air ↔ kelvin
block topology when the user does not name them explicitly, so .jou
files using a plain "concentric Kelvin" pattern need no manual sideset
work.

### BBND -- edge / curve labels (NGSolve `mesh.GetBBoundaries()`)

Source: Cubit **named curves** + Cubit **sidesets-on-curves**.

The exporter writes CD2 segments only for curves that own actual 1-D mesh edges
and whose parent surfaces both have exported Netgen descriptors.  Imported STL
or Sculpt geometry can retain CAD curves while producing no BBND mesh segment;
those stale curves are intentionally omitted from the companion edge-length
references.  For a point-only anchor, put the curve or vertices in a Cubit
**nodeset**; the exporter expands the nodeset to BBBND points (next table).

The BBND label-name convention to be respected once segment
generation lands:

| Cubit name on a curve | NGSolve BBND name | Meaning                                       |
|-----------------------|-------------------|-----------------------------------------------|
| `dir_<name>`          | `dir_<name>`      | Dirichlet edge (e.g. ground line in 2D)       |
| `neu_<name>`          | `neu_<name>`      | Explicit Neumann edge                         |
| anything else         | (passes through)  | Free-form; meaning is solver-side             |

### BBBND -- vertex / point labels (NGSolve `mesh.GetBBBoundaries()`)

Source: Cubit **nodesets**. Free-floating vertices (not merged into any
meshed volume, e.g. the bare vertex `add_kelvin_cubit` creates at the
Kelvin sphere centre) are anchored to the nearest mesh node so the
BBBND point is always usable as a Dirichlet anchor.

| Cubit nodeset name | NGSolve BBBND name | Meaning                                                    |
|--------------------|--------------------|------------------------------------------------------------|
| `GND`              | `GND`              | Special: Omega-reduced anchor at Kelvin sphere centre      |
| `dir_<name>`       | `dir_<name>`       | Dirichlet point (e.g. PEEC port gnd, source / sink reference) |
| anything else      | (passes through)   | Free-form name; meaning is solver-side                     |

`GND` is automatically created by the Auto-Kelvin helper at the
Kelvin sphere centre (the image of physical infinity) for use by
Omega-reduced FEM formulations.

## Python API

```python
import netgen          # must import before cubit (DLL load order)
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('open "model.cub5"')
cubit.cmd('mesh volume all')
cubit.cmd('block 1 add hex all')                          # elements must be in a block to export
cubit.cmd('export netgen "model.vol" order 3 overwrite')  # high-order CURVED .vol (order 1-5)

# Load it in NGSolve.  A high-order .vol already carries its curved mid-side nodes:
# load AS-IS and do NOT call mesh.Curve() -- mesh.Curve() re-curves from CAD geometry
# (absent in a loaded .vol) and would RESET every element to straight-sided.
from ngsolve import Mesh
mesh = Mesh("model.vol")
```

See `docs/cubit_mesh_export/hex_sphere_highorder/` for a runnable demo (a curved
hex sphere whose NGSolve volume converges to 4/3 pi r^3 as the order rises: order 1
-23 % -> order 2 -0.2 % -> order 3 +0.1 %).

The Cubit-side Python helpers (Kelvin transformation, etc.) live in
`cubit_mesh_export.cubit_helpers`:

```python
from cubit_mesh_export.cubit_helpers.add_kelvin import (
    add_kelvin_cubit,        # 3D Cubit path
    add_kelvin_occ,          # 3D OCC path
    add_kelvin_2d_axisym,    # 2D axisymmetric (r, z) path
    sym_sideset_name,        # canonical sym_<bc>_<axis> string
    parse_sym_label,         # inverse
)
```

In Cubit-embedded Python (where `cubit_mesh_export` itself is not
importable), the same helpers are available directly after
`cubit-plugin-install` deploys them to `<Cubit>/bin/plugins/cubit_helpers/`:

```python
# Inside a .jou or panel script, after add_kelvin is on sys.path:
python "import sys; sys.path.insert(0, r'<Cubit>/bin/plugins/cubit_helpers')"
python "from add_kelvin import add_kelvin_cubit"
python "add_kelvin_cubit(R=0.06, symmetry=['z'])"
```

The `export netgen ... add_kelvin` flow handles `sys.path`
itself, so users invoking Kelvin via the new APREPRO args do not need
to set anything by hand.

## Mesh consistency check (does NOT require Cubit)

```bash
check-vol model.vol                          # labels + curved-map quality
check-vol model.vol --strict-labels          # enforce canonical label names
check-vol model.vol --contract ih_labels.json --strict-labels
check-vol model.vol --json model.vol.json    # require this CAD reference
check-vol model.vol --tet-only --min-scaled-jacobian 0.05
check-vol model.vol --conductors copper,magnet \
  --sibc-boundaries conductor_air,conductor_exterior
check-vol model.vol --format json --report-json run/vol_check.json
```

The sibling `model.vol.json` written by `export netgen` is auto-discovered when
present. It is optional for a standalone `.vol`: mesh loading, label checks,
and the curved NGSolve mapping gate still run without Cubit or CAD data. Passing
`--json` makes that specific sidecar mandatory. Curved-map sampling is enabled
by default; `--no-quality` is available only for a quick label/CAD inspection.
`materials` contains only labels with a real CAD-volume reference;
`mesh_only_materials` contains free-mesh block labels that cannot supply one.
The checker rejects a label appearing in both sets, so a true CAD zero is not
confused with the absence of CAD ownership.

```python
from cubit_mesh_export.check import (
    check_consistency,
    check_label_contract,
    check_mesh_quality,
)

results = check_consistency(
    "model.vol",
    contract="ih_labels.json",
    strict_labels=True,
)
quality = check_mesh_quality(
    "model.vol",
    conductors=("copper", "magnet"),
    sibc_boundaries=("conductor_air", "conductor_exterior"),
    tet_only=True,
)
```

An application label contract is a versioned JSON object. `required` catches
missing solver labels; optional `allowed` lists reject labels that do not belong
to the selected application/mode.

```json
{
  "schema": "radia.vol-label-contract.v1",
  "application": "radia-ih/fem-kelvin",
  "strict_labels": true,
  "required": {
    "materials": ["coil", "air", "kelvin"],
    "boundaries": ["source", "sink", "sibc", "kelvin_int", "kelvin_ext"],
    "bbboundaries": ["GND"]
  },
  "allowed": {
    "materials": ["coil", "air", "kelvin"],
    "boundaries": [
      "source", "sink", "sibc", "coil_surface", "air_seam",
      "kelvin_int", "kelvin_ext"
    ]
  }
}
```

Strict mode requires descriptive lower snake-case material/boundary names,
while preserving Radia's `sym_bn=0_x` / `sym_ht=0_x` convention and the reserved
`GND` point label. It rejects generated fallbacks such as `volume_1` and
`Surface_7`, case-insensitive collisions, incomplete `source`/`sink` or
`kelvin_int`/`kelvin_ext` pairs, invalid Kelvin anchors, and contradictory
symmetry labels.

The quality gate samples the actual curved NGSolve element mapping; it does not
infer quality from straight corner nodes. A consistently positive or negative
element orientation is valid, while a sign change inside one element fails;
scaled quality uses `abs(det(J))`. It checks mapping magnitude, geometry order,
tetrahedron-only contracts, required labels, and material-aware face roles.
Only conductor-air or conductor-exterior faces may
be classified as SIBC.  Conductor-insulator faces retain a trace role, while
conductor-conductor faces retain the interface/loop-bridge role needed by the
reduced HCurl cycle space.

`check-vol` returns 0 for pass, 1 for a validation finding, and 2 for an input or
configuration error. JSON reports use schema
`cubit-mesh-export.vol-check.v1`. Run the checker after `.vol` export and before
solver or Simulink initialization. Material constants are not inferred from
mesh labels; the application's checked DesignSpec/configuration owns those
values and validates them separately.

## Part of the Radia project

Source: [github.com/ksugahar/Radia](https://github.com/ksugahar/Radia)
