# Mesh Format Routing: `.vol` vs `.step` — and why labels live in `.vol`

This document explains the two file formats that radia-* domain tools
accept as input, why labels live in `.vol`, and why
[`cubit-mesh-export`](../../packages/cubit-mesh-export/README.md) is
the single recommended chokepoint for producing labeled `.vol` files.

## TL;DR

- `.vol` (Netgen text format) carries **mesh + labels** at three
  dimension levels (material / BND / BBBND). Use it for any FEM
  workflow.
- `.step` (CAD geometry exchange) carries **geometry only**. Names
  on faces / solids are technically possible (AP242 + OCAF) but
  fragile in practice -- NGSolve OCC ignores them, most CAD tools
  drop them on round-trip, and there is no `sideset` / `nodeset`
  hierarchy. Use it only for label-free workflows (PEEC filament
  centerline extraction).
- `cubit-mesh-export` exists to be the **single chokepoint** for
  producing labeled `.vol`. STEP-with-extra-tricks alternatives are
  uniformly worse than just routing through Cubit.

## Format capability comparison

| Capability                                       | STEP (`.step`)                                   | Netgen Vol (`.vol`)                                  |
|--------------------------------------------------|--------------------------------------------------|------------------------------------------------------|
| Geometry (face / solid / edge topology)          | ✅ AP203 / AP214 / AP242 standard                 | (mesh derived from geometry; geometry detached after) |
| Single-face / single-solid name                  | △ AP242 + OCAF only; tool support inconsistent  | ✅ block / sideset / nodeset names                    |
| Group multiple faces under one label             | ❌ no representation                              | ✅ a sideset is by definition a face set               |
| Material / boundary / point hierarchy distinction | ❌ everything is "a name on a face/solid"         | ✅ four levels: material / BND / BBND / BBBND          |
| Symmetry / Kelvin / Dirichlet semantic conventions | ❌ would need a private string regex             | ✅ first-class: `sym_bn=0_x`, `kelvin_int`, `dir_<n>`  |
| High-order curving                                | ❌ linear only (CAD geometry, not mesh)          | ✅ p = 1..5 with ACIS projection                       |
| Tool portability                                  | △ names are fragile; many CAD tools drop them   | ✅ Netgen / NGSolve native, `Mesh("x.vol")` one-shot  |

## Operational consequences

### PEEC (`.step` is enough — label-free path)

PEEC works directly off `.step` because the workflow can derive
everything it needs from raw geometry topology:

- **Filament centerline** is extracted by tracing the longest
  `BSPLINE` chain through the coil solid (see
  `radia.coil_from_cad`).
- **Source / sink ports** can be identified geometrically
  (smallest-area faces at the wire ends).
- **No material / BC labels needed** -- the only physical input is
  the prescribed port current.

Driver:
```
python -m radia.radia_ih  --step coil.step      # PEEC mode
calc_peec.py --step coil.step
```

This works only because PEEC's input requirements are minimal.
The moment you need material distinction (multi-conductor
windings), surface impedance regions, or Kelvin open boundaries,
you are forced to leave `.step` and use `.vol`.

### FEM (`.vol` mandatory — labels carry the BC physics)

Every FEM panel reads `.vol` because it needs labels:

- Material name → assigns reluctivity / conductivity / source M.
- BND sideset name → assigns Dirichlet / Neumann / Robin / Kelvin /
  symmetry BC.
- BBBND nodeset name → assigns Omega-reduced anchor (e.g. `GND`
  vertex).

Driver:
```
python -m radia.radia_ih  model.vol             # FEM mode
python -m radia.radia_em  model.vol
calc_fem_kelvin.py --vol model.vol
calc_fem_coilmesh.py --vol model.vol
```

There is no `.step`-input FEM panel because there is no clean way
to attach the labels FEM needs to a STEP file.

### Hybrid (PEEC coil + FEM workpiece)

Mixed workflows like induction-heating PEEC+FEM (`calc_peec_bem.py`)
are intrinsically asymmetric:

- `--step coil.step`  → PEEC filament for the coil (label-free)
- `--vol  workpiece.vol` → FEM-SIBC for the workpiece (labels
  required: `workpiece` material, `sibc` boundary, optional
  `kelvin_int` / `kelvin_ext` for the open-boundary version)

This split is forced by what each format can carry, not by a
solver-side preference.

## Why `cubit-mesh-export` is the single chokepoint

Producing a labeled `.vol` requires a meshing tool that lets the
human (or an AI authoring `.jou`) attach **named blocks**, **named
sidesets**, and **named nodesets** to mesh entities. The lab's
production tool for this is Coreform Cubit, accessed through
[`cubit-mesh-export`](../../packages/cubit-mesh-export/README.md):

```
radia_export netgen "model.vol" order N \
    [add_kelvin] \
    [kelvin_mesh <m>] \
    [kelvin_sym_x|y|z {off|bn|ht}] \
    [overwrite]
```

Single command, single round trip: blocks → materials, sidesets →
BND, nodesets → BBBND, optional Kelvin geometry + sidesets, optional
symmetry-plane BC labels. Companion `.vol.json` is written next to
the `.vol` with CAD-reference Volume / Area / Length so a downstream
solver can sanity-check the mesh.

This is **the** reason cubit-mesh-export exists. Stripping it back
to "just the mesh export" loses the labeling layer that the rest of
the toolchain depends on.

## Why STEP-with-extra-tricks workarounds are not preferred

Three theoretical paths to "STEP + labels", all worse than routing
through `cubit-mesh-export`:

| Workaround                                           | Why it loses to `cubit-mesh-export`                                                                 |
|------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| OCC face name regex (`dir_gnd__face_3` etc.)         | CAD tools drop names on round-trip; humans + AIs find regex names hard to read; no group support    |
| STEP + sidecar JSON (`coil.step` + `coil.labels.json`) | Two-file pairing breaks if either is moved; no schema enforcement; reinvents `.vol` poorly         |
| STEP AP242 + OCAF document reader                    | NGSolve OCC ignores OCAF names; Netgen's reader skips them; partial Cubit support only             |

In every case the cost-benefit lands on the same answer: **if you
need labels, build the `.vol` from Cubit and stop trying**.

## Decision tree

```
Need labels (material / BND / BBBND) ?
├── No   → .step is enough.
│         Use it for PEEC, label-free CAD inspection, geometry sharing.
└── Yes  → .vol is required.
          Route through `cubit-mesh-export`:
              radia_export netgen "model.vol" order N [add_kelvin] [...]
          Then hand the .vol to the appropriate radia-* domain tool.
```

## See also

- [`cubit-mesh-export README`](../../packages/cubit-mesh-export/README.md)
  — APREPRO surface, label-name conventions (BND / BBND / BBBND tables),
  Kelvin / symmetry options.
- [`docs/cubit/export_NetgenMesh.md`](export_NetgenMesh.md) —
  details of the .vol writer.
- `radia_mcp.cubit.knowledge.format_routing` — same content as this
  document, exposed via `cubit_docs("format_routing")` for AI lookups.
