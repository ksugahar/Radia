"""Mesh format routing: .vol vs .step + why labels live in .vol.

Mirror of ``docs/cubit/Vol_vs_Step_Labels.md``.  Single source of
truth lives in the docs directory; this module is the AI-lookup
copy exposed via ``cubit_docs("format_routing")``.

Headline:
  - .vol carries mesh + labels (material / BND / BBBND).  FEM input.
  - .step carries geometry only.  PEEC input (label-free flow).
  - cubit-mesh-export is the single recommended chokepoint for
    producing labeled .vol; STEP-with-extra-tricks workarounds
    (OCC face name regex, sidecar JSON, AP242 OCAF reader) are
    uniformly worse.
"""

FORMAT_ROUTING = """
# Mesh Format Routing: .vol vs .step -- and why labels live in .vol

## TL;DR

- **`.vol`** (Netgen text) carries mesh + labels at three dimension
  levels (material / BND / BBBND).  Use it for any FEM workflow.
- **`.step`** (CAD geometry exchange) carries geometry only.  Names
  on faces / solids are technically possible (AP242 + OCAF) but
  fragile in practice -- NGSolve OCC ignores them, most CAD tools
  drop them on round-trip, and there is no `sideset` / `nodeset`
  hierarchy.  Use it only for label-free workflows (PEEC filament
  centerline extraction).
- **cubit-mesh-export** exists to be the single chokepoint for
  producing labeled `.vol`.  STEP-with-extra-tricks alternatives
  are uniformly worse than just routing through Cubit.

## Format capability comparison

| Capability                                  | STEP (`.step`)                                 | Netgen Vol (`.vol`)                                  |
|---------------------------------------------|------------------------------------------------|------------------------------------------------------|
| Geometry (face / solid / edge topology)     | YES (AP203 / AP214 / AP242)                   | (mesh derived from geometry; geometry detached after) |
| Single-face / single-solid name             | partial (AP242 + OCAF only; tool support inconsistent) | YES (block / sideset / nodeset names)         |
| Group multiple faces under one label        | NO (no representation)                         | YES (a sideset is by definition a face set)          |
| Material / boundary / point hierarchy       | NO (everything is "a name on a face/solid")    | YES (four levels: material / BND / BBND / BBBND)     |
| Symmetry / Kelvin / Dirichlet conventions   | NO (would need a private string regex)         | YES (first-class: sym_bn=0_x, kelvin_int, dir_<n>)   |
| High-order curving                          | NO (linear only -- it's CAD geometry, not mesh) | YES (p = 1..5 with ACIS projection)                 |
| Tool portability                            | partial (names fragile; many CAD tools drop them) | YES (Netgen / NGSolve native, Mesh("x.vol") one-shot) |

## Operational consequences

### PEEC: .step is enough (label-free path)

PEEC works directly off `.step` because the workflow can derive
everything it needs from raw geometry topology:

  - Filament centerline is extracted by tracing the longest BSPLINE
    chain through the coil solid (radia.coil_from_cad).
  - Source / sink ports can be identified geometrically (smallest-
    area faces at the wire ends).
  - No material / BC labels needed -- the only physical input is
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

### FEM: .vol mandatory (labels carry BC physics)

Every FEM panel reads `.vol` because it needs labels:

  - Material name -> reluctivity / conductivity / source M.
  - BND sideset name -> Dirichlet / Neumann / Robin / Kelvin /
    symmetry BC.
  - BBBND nodeset name -> Omega-reduced anchor (e.g. "GND" vertex).

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

Mixed workflows like induction-heating PEEC+FEM (calc_peec_bem.py)
are intrinsically asymmetric:

  --step coil.step      -> PEEC filament for the coil (label-free)
  --vol  workpiece.vol  -> FEM-SIBC for the workpiece (labels:
                            "workpiece" material, "sibc" boundary,
                            optional "kelvin_int" / "kelvin_ext"
                            for the open-boundary version)

This split is forced by what each format can carry, not by a
solver-side preference.

## Why cubit-mesh-export is the single chokepoint

Producing a labeled `.vol` requires a meshing tool that lets the
human (or an AI authoring `.jou`) attach **named blocks**, **named
sidesets**, and **named nodesets** to mesh entities.  The lab's
production tool for this is Coreform Cubit, accessed through
cubit-mesh-export:

```
export netgen "model.vol" order N \\
    [add_kelvin] \\
    [kelvin_mesh <m>] \\
    [kelvin_sym_x|y|z {off|bn|ht}] \\
    [overwrite]
```

Single command, single round trip:

  - blocks -> materials
  - sidesets -> BND
  - nodesets -> BBBND
  - optional Kelvin geometry + sidesets
  - optional symmetry-plane BC labels
  - companion `.vol.json` written next to the `.vol` with
    CAD-reference Volume / Area / Length so a downstream solver
    can sanity-check the mesh.

This is the reason cubit-mesh-export exists.  Stripping it back
to "just the mesh export" loses the labeling layer that the rest
of the toolchain depends on.

## Why STEP-with-extra-tricks workarounds are not preferred

Three theoretical paths to "STEP + labels", all worse than routing
through cubit-mesh-export:

| Workaround                                            | Why it loses to cubit-mesh-export                                                                  |
|-------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| OCC face name regex (`dir_gnd__face_3` etc.)          | CAD tools drop names on round-trip; humans + AIs find regex names hard to read; no group support   |
| STEP + sidecar JSON (`coil.step` + `coil.labels.json`) | Two-file pairing breaks if either is moved; no schema enforcement; reinvents `.vol` poorly       |
| STEP AP242 + OCAF document reader                     | NGSolve OCC ignores OCAF names; Netgen's reader skips them; partial Cubit support only             |

In every case the cost-benefit lands on the same answer: **if you
need labels, build the `.vol` from Cubit and stop trying**.

## Decision tree

```
Need labels (material / BND / BBBND) ?
+-- No   -> .step is enough.
|           Use it for PEEC, label-free CAD inspection, geometry sharing.
+-- Yes  -> .vol is required.
            Route through cubit-mesh-export:
                export netgen "model.vol" order N [add_kelvin] [...]
            Then hand the .vol to the appropriate radia-* domain tool.
```

## Mesh-type routing inside `.vol`

The `.vol` extension is a container, not a promise that the mesh is tri/tet.
Route by semantic inventory:

- **tet-only generation**: use build123d/Netgen/OCC for the first-order
  H1/HCurl/FEM-BEM education path.  The radia-ngsolve tri/tet parser rejects
  quad/hex/wedge/pyramid instead of silently converting them.
- **hex-led or mixed hex+pyramid+tet generation**: use Cubit/Coreform.  This is
  where Cubit adds unique value.  Run `cubit_vol_inventory` first; it reports
  triangle/quad surface records and tet/pyramid/wedge/hex volume records, then
  tells the agent whether the file is a `netgen_tri_tet_path` or a
  `cubit_hex_or_mixed_path`.
- **pyramid transition elements** are expected at some hex/tet interfaces.
  Inventory them explicitly and route to a solver/export contract that supports
  them, or perform a deliberate conversion.  Do not split them implicitly in a
  generic `.vol` reader.
- **high-order `.vol` files** still route by the first-order element arity in
  `surfaceelements` and `volumeelements`.  Order-2..5 exports add
  `curvedelements` data and grow the file, but the hex/pyramid/tet inventory
  remains the same routing contract.  Do not infer element existence from the
  companion `.vol.json` material volume alone: transition blocks such as
  pyramids can legitimately report zero material volume in the sidecar while
  still being present in the `.vol` topology.

## See also

- `docs/cubit/Vol_vs_Step_Labels.md` -- canonical (this file is a
  mirror for AI lookup).
- `packages/cubit-mesh-export/README.md` -- APREPRO surface,
  label-name conventions (BND / BBND / BBBND tables), Kelvin /
  symmetry options.
- `docs/cubit/export_NetgenMesh.md` -- details of the .vol writer.
"""


def get_format_routing_documentation(topic: str = "format_routing") -> str:
    """Return mesh-format-routing documentation by topic.

    Topics:
        "format_routing"   -- (default) the .vol-vs-.step routing rule:
                              .vol carries labels (FEM input), .step is
                              geometry-only (PEEC input), cubit-mesh-
                              export is the single labeled-.vol chokepoint.
    """
    topic = (topic or "").lower().strip()
    if topic in (
        "",
        "format_routing", "routing", "vol_vs_step",
        "format", "step_vs_vol", "labels_format",
    ):
        return FORMAT_ROUTING
    return (
        f"Unknown format-routing topic '{topic}'.  Available: "
        f"format_routing."
    )
